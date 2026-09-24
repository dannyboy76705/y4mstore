/*
 * y4mstore - store raw YUV frame data from a y4m stream and serve it
 *            back out as y4m; can also sort stored frames by nilsimsa
 *            similarity before archiving
 *
 * MIT License
 *
 * Copyright (c) 2026 Daniel Lee Witzel (nilsort)
 * Copyright (c) 2026 Daniel Lee Witzel (y4mstore)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * y4mstore is a spinoff of nilsort, narrowed to raw (uncompressed) YUV
 * frame files. The MinHash hasher and the EBML/Matroska container
 * handling (-M / -m) were removed: on raw YUV, plain nilsimsa measured
 * far better than MinHash in nilsort's benchmarks, and there is no
 * container framing to strip. y4mstore/y4m (yuv4mpeg) support is planned.
 *
 * The nilsimsa hashing algorithm implemented below (TRAN/POPC-derived
 * tables, tran3 mixing function, digest construction) is a faithful
 * port of the long-standing public reference implementation originally
 * written by cmeclax, based on Damiani et al. 2004, "An Open
 * Digest-based Technique for Spam Detection". The threading,
 * greedy + 2-opt ordering, progress reporting, and POPCNT-based
 * comparison are inherited from nilsort.
 *
 * Directory layout: a set directory holds specs.txt -- in y4m terms:
 * resolution=640x480 / sampling=420 (or 422, 444) / framerate=24000:1001 /
 * bitdepth=8 (or 10) / aspect=1:1 -- describing the raw YUV common to
 * every clip, plus one subdirectory per clip. A directory that holds
 * frame files directly is a clip: its specs come from the parent (set)
 * directory, whichever mode it is named in. If specs.txt is missing,
 * y4mstore asks for the first four values on the terminal and writes it;
 * aspect is never asked for (footage is resampled to square pixels first)
 * and is added automatically, with a message, if the line is absent. (Reading a file list from
 * stdin with "-" has no directory, so no specs.txt is involved.)
 *
 * Serve mode (-s): "y4mstore -s <clipdir>" streams a clip directory (or a
 * directory of y/u/v planes written by -e, the reverse of -e) as
 * yuv4mpeg (y4m) on stdout -- or to -f FILE -- so it can be piped to
 * ffplay, mpv or ffmpeg, or saved as a .y4m file. The stream header is
 * built from specs.txt in the clip's parent directory; the raw frame
 * files are already in y4m plane order (Y, U, V) and are copied through.
 *
 * Several directories can be sorted together (sorting only): all of them
 * whole-frame clips, or all of them y/u/v component directories (checked
 * first). One list comes out: for frames, every frame file pooled and sorted
 * once; for components, all the y planes sorted, then all the u, then all the
 * v. Directories are taken in the order given.
 *
 * A directory that already holds y/, u/ and v/ planes (what -e writes) is
 * recognised without -e: the selected sorter (nilsimsa, or -a) runs on each
 * component directory and the three lists are appended, y then u then v.
 *
 * Scans a directory (or reads a file list from stdin), computes a
 * nilsimsa locality-sensitive hash digest for each file, then
 * greedily orders the files so that similar files end up adjacent
 * in the output list. The ordered file list goes to stdout, so you can
 * pipe it into an archiver the same way you would with binsort:
 *
 *   $ y4mstore -n <dir> | tar -T- --no-recursion -czf out.tar.gz
 *
 * With `-n`, this tool treats its ENTIRE input as one batch, so every
 * file is weighed against every other file -- e.g. point it at one
 * episode's frame directory and all of that episode's frames are
 * compared against each other in a single sort. This maximizes accuracy
 * at the cost of O(n^2) time. Memory is not the constraint here: digests
 * are only 32 bytes each, so even a million files is ~32 MB of digest
 * storage. If you have multiple episodes, run y4mstore once per episode
 * directory rather than once across all episodes combined, unless you
 * deliberately want the pooled cross-episode sort described below.
 *
 * On top of the initial greedy nearest-neighbor construction, this
 * tool also runs 2-opt local-search refinement passes by default
 * (like binsort's -o optimization level, but here defaulting to
 * "run until it stops improving" since accuracy matters more than
 * speed). Each pass is another O(n^2) scan; use -o to cap or disable
 * it if a run is taking too long.
 *
 * Build:
 *   gcc -O2 -Wall -pthread -march=native -o y4mstore y4mstore.c
 * (-march=native builds for the exact CPU you're compiling on and
 * picks up hardware POPCNT for digest comparison, ~8x faster than a
 * lookup-table popcount in testing. To build once and run on other
 * x86-64 machines from ~2008 (Intel) / ~2012 (AMD) or later, use
 * -mpopcnt instead. Without either flag gcc emits a portable software
 * popcount: correct everywhere, just slower.)
 *
 * Usage:
 *   y4mstore [options] <dir>
 *   find /some/dir -type f | y4mstore [options] -
 *
 * Without -n or -a, no sort runs at all: the file list is printed in
 * natural filename order (no hashing, no O(n^2) comparisons, no
 * duplicate detection -- that's paired with -n, see below). This is
 * the default so that adding clips over time, or just checking a
 * growing pooled set, doesn't force a full nilsimsa sort on every run;
 * pass -n once, when you're ready to actually sort and archive.
 *
 * Options:
 *   -n     nilsimsa similarity sort (this used to be the unconditional
 *          default). Hashes every file, greedily chains them by
 *          similarity, then refines with 2-opt. Exact-duplicate
 *          detection (a byte-identical-file check) is a byproduct of
 *          this hashing pass, so it only runs with -n, not by default.
 *   -o N   max nilsimsa passes (default 50; -n only). 0 disables the
 *          2-opt refinement stage.
 *   -t N   number of hashing threads (-n), or averaging threads (-a),
 *          default = number of CPUs
 *   -f PATH write the file list to PATH instead of stdout
 *          (with -s: write the y4m stream to PATH)
 *   -s     serve the named clip directory as a y4m stream (see above)
 *   -e SRC split a y4m stream into y/u/v plane files (one directory each), then
 *          sort each plane directory (with -n or -a; default: natural order,
 *          no sort) into one list
 *   -i SRC ingest a y4m stream (file, or - for stdin) into a new clip directory of
 *          raw frame files, creating/checking specs.txt from its header
 *   -a     alternate sort: order frames by their additive Y/U/V component
 *          averages (scaled x100 for precision; a bit-depth-dependent
 *          digit width per plane) instead of nilsimsa similarity.
 *          Mutually exclusive with -n.
 *   -c     check the set directory's specs.txt: create it if missing
 *          (asking for the values), otherwise list it; nothing else runs
 *   -d PATH with -a: outputs averaged per frame/component to PATH
 *   -p     with -i/-e reading stdin: a running count of frames read (opt-in)
 *   -q     quiet (suppress progress messages on stderr)
 *   -h     show this help
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <ctype.h>
#include <strings.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <stdatomic.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <stdarg.h>

/* ---------------------------------------------------------------- */
/* Progress bar (stderr only, throttled to ~10 updates/sec)         */
/* ---------------------------------------------------------------- */

static void progress_bar(const char *label, size_t done, size_t total) {
    static struct timespec last_print;
    static int have_last = 0;
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    if (have_last && done != total) {
        double dt = (double)(now.tv_sec - last_print.tv_sec) +
                    (double)(now.tv_nsec - last_print.tv_nsec) / 1e9;
        if (dt < 0.1) return; /* throttle so huge n doesn't flood stderr */
    }
    last_print = now;
    have_last = 1;

    const int width = 30;
    double frac = total ? (double)done / (double)total : 1.0;
    if (frac > 1.0) frac = 1.0;
    int filled = (int)(frac * width);

    fprintf(stderr, "\r%-10s [", label);
    for (int i = 0; i < width; i++) fputc(i < filled ? '#' : '-', stderr);
    fprintf(stderr, "] %5.1f%% (%zu/%zu)", frac * 100.0, done, total);
    if (done >= total) fputc('\n', stderr);
    fflush(stderr);
}

/* ---------------------------------------------------------------- */
/* Nilsimsa core                                                    */
/* ---------------------------------------------------------------- */

static const unsigned char TRAN[256] = {
    0x02,0xD6,0x9E,0x6F,0xF9,0x1D,0x04,0xAB,0xD0,0x22,0x16,0x1F,0xD8,0x73,0xA1,0xAC,
    0x3B,0x70,0x62,0x96,0x1E,0x6E,0x8F,0x39,0x9D,0x05,0x14,0x4A,0xA6,0xBE,0xAE,0x0E,
    0xCF,0xB9,0x9C,0x9A,0xC7,0x68,0x13,0xE1,0x2D,0xA4,0xEB,0x51,0x8D,0x64,0x6B,0x50,
    0x23,0x80,0x03,0x41,0xEC,0xBB,0x71,0xCC,0x7A,0x86,0x7F,0x98,0xF2,0x36,0x5E,0xEE,
    0x8E,0xCE,0x4F,0xB8,0x32,0xB6,0x5F,0x59,0xDC,0x1B,0x31,0x4C,0x7B,0xF0,0x63,0x01,
    0x6C,0xBA,0x07,0xE8,0x12,0x77,0x49,0x3C,0xDA,0x46,0xFE,0x2F,0x79,0x1C,0x9B,0x30,
    0xE3,0x00,0x06,0x7E,0x2E,0x0F,0x38,0x33,0x21,0xAD,0xA5,0x54,0xCA,0xA7,0x29,0xFC,
    0x5A,0x47,0x69,0x7D,0xC5,0x95,0xB5,0xF4,0x0B,0x90,0xA3,0x81,0x6D,0x25,0x55,0x35,
    0xF5,0x75,0x74,0x0A,0x26,0xBF,0x19,0x5C,0x1A,0xC6,0xFF,0x99,0x5D,0x84,0xAA,0x66,
    0x3E,0xAF,0x78,0xB3,0x20,0x43,0xC1,0xED,0x24,0xEA,0xE6,0x3F,0x18,0xF3,0xA0,0x42,
    0x57,0x08,0x53,0x60,0xC3,0xC0,0x83,0x40,0x82,0xD7,0x09,0xBD,0x44,0x2A,0x67,0xA8,
    0x93,0xE0,0xC2,0x56,0x9F,0xD9,0xDD,0x85,0x15,0xB4,0x8A,0x27,0x28,0x92,0x76,0xDE,
    0xEF,0xF8,0xB2,0xB7,0xC9,0x3D,0x45,0x94,0x4B,0x11,0x0D,0x65,0xD5,0x34,0x8B,0x91,
    0x0C,0xFA,0x87,0xE9,0x7C,0x5B,0xB1,0x4D,0xE5,0xD4,0xCB,0x10,0xA2,0x17,0x89,0xBC,
    0xDB,0xB0,0xE2,0x97,0x88,0x52,0xF7,0x48,0xD3,0x61,0x2C,0x3A,0x2B,0xD1,0x8C,0xFB,
    0xF1,0xCD,0xE4,0x6A,0xE7,0xA9,0xFD,0xC4,0x37,0xC8,0xD2,0xF6,0xDF,0x58,0x72,0x4E
};

/* A nilsimsa digest is 256 bits. */
#define DIGEST_BYTES 32

/* Similarity score via hardware POPCNT: compare the two 32-byte digests
 * as four 64-bit words, XOR, and count set bits with __builtin_popcountll
 * (compiles to a single POPCNT instruction with -mpopcnt / -march=native
 * on x86, or a portable software fallback otherwise -- correct either
 * way, just faster when the hardware instruction is available). This
 * replaces a per-byte lookup-table approach that measured ~8x slower
 * for this same computation on typical x86 hardware; note this affects
 * only *how fast* the comparison runs, not the digests or scores
 * themselves, so output is bit-for-bit identical either way.
 * 128 = identical digests, going down as they diverge. */
static inline int nilsimsa_compare(const unsigned char *d1, const unsigned char *d2) {
    uint64_t a[4], b[4];
    memcpy(a, d1, 32);
    memcpy(b, d2, 32);
    int bits = __builtin_popcountll(a[0] ^ b[0]) +
               __builtin_popcountll(a[1] ^ b[1]) +
               __builtin_popcountll(a[2] ^ b[2]) +
               __builtin_popcountll(a[3] ^ b[3]);
    return 128 - bits;
}

typedef struct {
    long count;         /* bytes seen so far */
    uint32_t acc[256];  /* trigram accumulators */
    int lastch[4];      /* last 4 bytes seen, -1 = not yet valid */
} nilsimsa_t;

static void nilsimsa_init(nilsimsa_t *n) {
    n->count = 0;
    memset(n->acc, 0, sizeof(n->acc));
    n->lastch[0] = n->lastch[1] = n->lastch[2] = n->lastch[3] = -1;
}

static inline unsigned char tran3(unsigned char a, unsigned char b, unsigned char c, int nn) {
    return (unsigned char)(((TRAN[(a + nn) & 255] ^ (TRAN[b] * (nn + nn + 1))) +
                             TRAN[c ^ TRAN[nn]]) & 255);
}

static void nilsimsa_update(nilsimsa_t *n, const unsigned char *data, size_t len) {
    int *l = n->lastch;
    for (size_t k = 0; k < len; k++) {
        unsigned char ch = data[k];
        n->count++;
        if (l[1] > -1) {
            n->acc[tran3(ch, (unsigned char)l[0], (unsigned char)l[1], 0)]++;
        }
        if (l[2] > -1) {
            n->acc[tran3(ch, (unsigned char)l[0], (unsigned char)l[2], 1)]++;
            n->acc[tran3(ch, (unsigned char)l[1], (unsigned char)l[2], 2)]++;
        }
        if (l[3] > -1) {
            n->acc[tran3(ch, (unsigned char)l[0], (unsigned char)l[3], 3)]++;
            n->acc[tran3(ch, (unsigned char)l[1], (unsigned char)l[3], 4)]++;
            n->acc[tran3(ch, (unsigned char)l[2], (unsigned char)l[3], 5)]++;
            n->acc[tran3((unsigned char)l[3], (unsigned char)l[0], ch, 6)]++;
            n->acc[tran3((unsigned char)l[3], (unsigned char)l[2], ch, 7)]++;
        }
        l[3] = l[2]; l[2] = l[1]; l[1] = l[0]; l[0] = ch;
    }
}

/* Produce the 32-byte (256-bit) digest from accumulated state. */
static void nilsimsa_digest(const nilsimsa_t *n, unsigned char out[32]) {
    long total;
    if (n->count == 3)      total = 1;
    else if (n->count == 4) total = 4;
    else if (n->count > 4)  total = 8 * n->count - 28;
    else                    total = 0;

    double threshold = total / 256.0;
    unsigned char code[32];
    memset(code, 0, sizeof(code));
    for (int i = 0; i < 256; i++) {
        if ((double)n->acc[i] > threshold) {
            code[i >> 3] |= (unsigned char)(1u << (i & 7));
        }
    }
    /* reference implementation reverses byte order in the result */
    for (int i = 0; i < 32; i++) out[i] = code[31 - i];
}

/* ---------------------------------------------------------------- */
/* Exact-duplicate fingerprint                                       */
/* ---------------------------------------------------------------- */

/* 64-bit FNV-1a, computed alongside the nilsimsa digest over exactly
 * the same bytes (the whole file). This is for EXACT duplicate detection -- nilsimsa is
 * a fuzzy/locality-sensitive hash and isn't meant to prove equality,
 * only similarity. Paired with a length check, FNV-1a collisions on
 * real (non-adversarial) file content are negligible in practice. */
#define FNV_OFFSET_64 0xcbf29ce484222325ULL
#define FNV_PRIME_64  0x100000001b3ULL

static inline uint64_t fnv1a64_update(uint64_t h, const unsigned char *data, size_t len) {
    for (size_t i = 0; i < len; i++) {
        h ^= data[i];
        h *= FNV_PRIME_64;
    }
    return h;
}

static int hash_file(const char *path, unsigned char *digest_out,
                      uint64_t *dup_hash_out, uint64_t *dup_len_out) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    nilsimsa_t n;
    nilsimsa_init(&n);
    uint64_t fh = FNV_OFFSET_64;
    uint64_t total = 0;
    unsigned char buf[1 << 16];
    size_t r;
    while ((r = fread(buf, 1, sizeof(buf), f)) > 0) {
        nilsimsa_update(&n, buf, r);
        fh = fnv1a64_update(fh, buf, r);
        total += r;
    }
    fclose(f);
    nilsimsa_digest(&n, digest_out);
    *dup_hash_out = fh;
    *dup_len_out = total;
    return 0;
}

typedef struct {
    char **paths;
    unsigned char *digests;
    uint64_t *dup_hash;
    uint64_t *dup_len;
    size_t start, end;
    int quiet;
    atomic_size_t *done_counter;
} hash_job_t;

static void *hash_worker(void *arg) {
    hash_job_t *j = (hash_job_t *)arg;
    for (size_t i = j->start; i < j->end; i++) {
        unsigned char *slot = j->digests + i * DIGEST_BYTES;
        int rc = hash_file(j->paths[i], slot, &j->dup_hash[i], &j->dup_len[i]);
        if (rc != 0) {
            if (!j->quiet) {
                fprintf(stderr, "\ny4mstore: warning: cannot read %s: %s\n",
                        j->paths[i], strerror(errno));
            }
            memset(slot, 0, DIGEST_BYTES);
            j->dup_hash[i] = 0;
            j->dup_len[i] = 0;
        }
        atomic_fetch_add(j->done_counter, 1);
    }
    return NULL;
}

static void hash_all(char **paths, size_t n, unsigned char *digests,
                      uint64_t *dup_hash, uint64_t *dup_len,
                      int nthreads, int quiet) {
    if (n == 0) return;
    if (nthreads < 1) nthreads = 1;
    if ((size_t)nthreads > n) nthreads = (int)n;

    atomic_size_t done_counter = 0;
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)nthreads);
    hash_job_t *jobs = malloc(sizeof(hash_job_t) * (size_t)nthreads);
    size_t chunk = (n + (size_t)nthreads - 1) / (size_t)nthreads;
    int launched = 0;

    for (int t = 0; t < nthreads; t++) {
        size_t start = (size_t)t * chunk;
        if (start >= n) break;
        size_t end = start + chunk;
        if (end > n) end = n;
        jobs[t].paths = paths;
        jobs[t].digests = digests;
        jobs[t].dup_hash = dup_hash;
        jobs[t].dup_len = dup_len;
        jobs[t].start = start;
        jobs[t].end = end;
        jobs[t].quiet = quiet;
        jobs[t].done_counter = &done_counter;
        pthread_create(&th[t], NULL, hash_worker, &jobs[t]);
        launched++;
    }

    if (!quiet) {
        size_t seen;
        do {
            struct timespec nap = {0, 50 * 1000 * 1000}; /* 50ms */
            nanosleep(&nap, NULL);
            seen = atomic_load(&done_counter);
            progress_bar("hashing", seen, n);
        } while (seen < n);
    }

    for (int t = 0; t < launched; t++) pthread_join(th[t], NULL);
    free(th);
    free(jobs);
}

/* ---------------------------------------------------------------- */
/* Exact duplicate detection (over the same bytes hashed above)     */
/* ---------------------------------------------------------------- */

typedef struct {
    size_t idx;
    uint64_t len;
    uint64_t fh;
} dup_key_t;

static int dup_key_cmp(const void *a, const void *b) {
    const dup_key_t *ka = a, *kb = b;
    if (ka->len != kb->len) return (ka->len < kb->len) ? -1 : 1;
    if (ka->fh != kb->fh) return (ka->fh < kb->fh) ? -1 : 1;
    return 0;
}

/* Groups files by (content length, FNV-1a fingerprint). Prints a
 * one-line summary (respecting quiet) of how many files are exact
 * duplicates of at least one other file, and how many duplicate sets
 * that forms. Does not print individual file paths or group contents
 * -- just the totals the person asked for. */
static void report_duplicates(uint64_t *dup_hash, uint64_t *dup_len,
                               size_t n, int quiet) {
    if (n < 2) {
        if (!quiet) fprintf(stderr, "y4mstore: 0 duplicate files (0 sets)\n");
        return;
    }
    dup_key_t *keys = malloc(sizeof(dup_key_t) * n);
    for (size_t i = 0; i < n; i++) {
        keys[i].idx = i;
        keys[i].len = dup_len[i];
        keys[i].fh = dup_hash[i];
    }
    qsort(keys, n, sizeof(dup_key_t), dup_key_cmp);

    size_t dup_files = 0, dup_sets = 0;
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && keys[j].len == keys[i].len && keys[j].fh == keys[i].fh) j++;
        size_t group_size = j - i;
        if (group_size > 1) {
            dup_files += group_size;
            dup_sets++;
        }
        i = j;
    }
    free(keys);

    if (!quiet) {
        fprintf(stderr, "y4mstore: %zu duplicate file%s found in %zu set%s of identical content\n",
                dup_files, dup_files == 1 ? "" : "s",
                dup_sets, dup_sets == 1 ? "" : "s");
    }
}

/* ---------------------------------------------------------------- */
/* Greedy similarity ordering (O(n^2) comparisons, O(n) memory)     */
/* ---------------------------------------------------------------- */

static inline const unsigned char *digest_at(const unsigned char *digests, size_t i) {
    return digests + i * DIGEST_BYTES;
}

static size_t *greedy_order(unsigned char *digests, size_t n, int quiet) {
    if (n == 0) return NULL;
    size_t *order = malloc(sizeof(size_t) * n);
    unsigned char *visited = calloc(n, 1);

    order[0] = 0;
    visited[0] = 1;
    size_t cur = 0;

    for (size_t step = 1; step < n; step++) {
        int best_score = -1000;
        size_t best_idx = 0;
        int found = 0;
        for (size_t i = 0; i < n; i++) {
            if (visited[i]) continue;
            int s = nilsimsa_compare(digest_at(digests, cur),
                                    digest_at(digests, i));
            if (!found || s > best_score) {
                best_score = s;
                best_idx = i;
                found = 1;
            }
        }
        order[step] = best_idx;
        visited[best_idx] = 1;
        cur = best_idx;
        if (!quiet) progress_bar("ordering", step, n - 1);
    }
    free(visited);
    return order;
}

/* One full 2-opt sweep over the current order. For every pair of
 * non-adjacent "edges" (order[i],order[i+1]) and (order[j],order[j+1]),
 * check whether reversing the segment between them increases the
 * total similarity of the two boundary edges; if so, apply it
 * immediately (first-improvement strategy) and keep scanning.
 * Returns 1 if any improving move was applied this sweep, else 0.
 */
static int twoopt_pass(unsigned char *digests, size_t *order, size_t n,
                        const char *label, int quiet) {
    int improved = 0;
    if (n < 4) return 0;
    size_t last_i = n - 4; /* i ranges 0 .. n-4 inclusive (i+3 < n) */
    for (size_t i = 0; i + 3 < n; i++) {
        int score_i = nilsimsa_compare(digest_at(digests, order[i]),
                                      digest_at(digests, order[i + 1]));
        for (size_t j = i + 2; j + 1 < n; j++) {
            int old_score = score_i +
                nilsimsa_compare(digest_at(digests, order[j]),
                                digest_at(digests, order[j + 1]));
            int new_score =
                nilsimsa_compare(digest_at(digests, order[i]),
                                digest_at(digests, order[j])) +
                nilsimsa_compare(digest_at(digests, order[i + 1]),
                                digest_at(digests, order[j + 1]));
            if (new_score > old_score) {
                size_t lo = i + 1, hi = j;
                while (lo < hi) {
                    size_t tmp = order[lo];
                    order[lo] = order[hi];
                    order[hi] = tmp;
                    lo++; hi--;
                }
                improved = 1;
                score_i = nilsimsa_compare(digest_at(digests, order[i]),
                                          digest_at(digests, order[i + 1]));
            }
        }
        if (!quiet) progress_bar(label, i, last_i);
    }
    return improved;
}

/* Run up to max_sweeps full 2-opt sweeps, stopping early once a sweep
 * makes no improving move (i.e. it has converged to a local optimum). */
static void optimize_order(unsigned char *digests, size_t *order,
                            size_t n, int max_sweeps, int quiet) {
    if (max_sweeps <= 0) return;
    char label[24];
    for (int s = 0; s < max_sweeps; s++) {
        snprintf(label, sizeof(label), "sweep %d", s + 1);
        int improved = twoopt_pass(digests, order, n, label, quiet);
        if (!improved) {
            if (!quiet) {
                fprintf(stderr, "y4mstore: converged after %d sweep(s)\n", s + 1);
            }
            break;
        }
    }
}

/* ---------------------------------------------------------------- */
/* specs.txt: raw YUV properties shared by every clip in a set       */
/*                                                                   */
/* Directory layout:                                                 */
/*     <set>/specs.txt      short file, common to every clip         */
/*     <set>/<clip>/...     one subdirectory per clip                */
/*                                                                   */
/* specs.txt is key=value lines, in yuv4mpeg (y4m) terms:            */
/*     resolution=640x480      -> W640 H480                          */
/*     sampling=420            -> C420jpeg  (C422, C444)             */
/*     framerate=24000:1001    -> F24000:1001                        */
/*     bitdepth=8              -> no suffix (10 -> C420p10 etc)      */
/*     aspect=1:1              -> A1:1  (never asked for; see below) */
/* Blank lines and lines starting with '#' are ignored. If it is     */
/* missing, the user is asked for the first four values and it is    */
/* written. aspect (pixel aspect ratio) is always 1:1 -- footage is  */
/* resampled to square pixels before storing -- so it is never asked */
/* for; a specs.txt without an aspect line gets "aspect=1:1"         */
/* appended, with a message.                                         */
/* ---------------------------------------------------------------- */

#define SPECS_FILENAME "specs.txt"

typedef struct {
    const char *name;    /* value used in specs.txt */
    const char *desc;    /* shown in the prompt */
    int hsub, vsub;      /* chroma subsampling divisors (horizontal, vertical) */
    const char *y4m8;    /* y4m C tag at 8 bits */
    const char *y4m10;   /* y4m C tag at 10 bits */
} sampling_info_t;

/* Planar Y, U, V -- the plane order y4m always uses. */
static const sampling_info_t SAMPLINGS[] = {
    {"420", "4:2:0 - chroma half width, half height", 2, 2, "420jpeg", "420p10"},
    {"422", "4:2:2 - chroma half width, full height", 2, 1, "422",     "422p10"},
    {"444", "4:4:4 - chroma full width, full height", 1, 1, "444",     "444p10"},
};
#define N_SAMPLINGS (sizeof(SAMPLINGS) / sizeof(SAMPLINGS[0]))

/* Bits per sample. 8 is one byte per sample; 10 is a 16-bit little-endian
 * word per sample, as in y4m's p10 colorspaces (yuv420p10le). */
static const int BITDEPTHS[] = {8, 10};
#define N_BITDEPTHS (sizeof(BITDEPTHS) / sizeof(BITDEPTHS[0]))

typedef struct {
    int width, height;
    int sampling;            /* index into SAMPLINGS */
    long fps_num, fps_den;   /* frame rate as a reduced fraction */
    int bitdepth;            /* one of BITDEPTHS */
    int have_aspect;         /* specs.txt already had an aspect line */
    uint64_t frame_bytes;    /* bytes in one raw frame at these specs */
} specs_t;

static char *trim(char *s) {
    while (isspace((unsigned char)*s)) s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1])) *--e = 0;
    return s;
}

/* Lets an answer be typed either as the bare value ("640x480") or as the
 * whole example line ("resolution=640x480"). */
static char *strip_key_prefix(char *s, const char *key) {
    size_t kl = strlen(key);
    if (strncasecmp(s, key, kl) == 0) {
        char *p = s + kl;
        while (isspace((unsigned char)*p)) p++;
        if (*p == '=') {
            p++;
            while (isspace((unsigned char)*p)) p++;
            return p;
        }
    }
    return s;
}

/* "640x480" (spaces around the x are tolerated). Returns 0 on success. */
static int parse_resolution(const char *s, int *w, int *h) {
    char *end;
    long a = strtol(s, &end, 10);
    if (end == s) return -1;
    while (isspace((unsigned char)*end)) end++;
    if (*end != 'x' && *end != 'X') return -1;
    const char *p = end + 1;
    long b = strtol(p, &end, 10);
    if (end == p) return -1;
    while (isspace((unsigned char)*end)) end++;
    if (*end != 0) return -1;
    if (a < 1 || a > 65535 || b < 1 || b > 65535) return -1;
    *w = (int)a;
    *h = (int)b;
    return 0;
}

/* Returns an index into SAMPLINGS, or -1. Exactly "420", "422" or "444". */
static int parse_sampling(const char *s) {
    for (size_t i = 0; i < N_SAMPLINGS; i++) {
        if (!strcmp(s, SAMPLINGS[i].name)) return (int)i;
    }
    return -1;
}

static long gcd_l(long a, long b) {
    while (b) { long t = a % b; a = b; b = t; }
    return a;
}

/* Fraction form only, in y4m's notation "24000:1001" (a "/" is accepted
 * too: "24000/1001"). Decimals ("29.97") and bare integers ("24") are
 * rejected on purpose -- an exact N:D is unambiguous. Stored reduced
 * (48000:2002 -> 24000:1001). Returns 0 on success. */
static int parse_fps(const char *s, long *num, long *den) {
    char *end;
    long n = strtol(s, &end, 10);
    if (end == s) return -1;
    while (isspace((unsigned char)*end)) end++;
    if (*end != ':' && *end != '/') return -1;
    const char *p = end + 1;
    long d = strtol(p, &end, 10);
    if (end == p) return -1;
    while (isspace((unsigned char)*end)) end++;
    if (*end != 0) return -1;
    if (n < 1 || d < 1) return -1;
    if ((double)n / (double)d > 1000.0) return -1;
    long g = gcd_l(n, d);
    *num = n / g;
    *den = d / g;
    return 0;
}

/* 8 or 10. Returns 0 on success. */
static int parse_bitdepth(const char *s, int *out) {
    char *end;
    long v = strtol(s, &end, 10);
    if (end == s) return -1;
    while (isspace((unsigned char)*end)) end++;
    if (*end != 0) return -1;
    for (size_t i = 0; i < N_BITDEPTHS; i++) {
        if (v == BITDEPTHS[i]) { *out = (int)v; return 0; }
    }
    return -1;
}

/* Pixel aspect ratio: only square pixels, "1:1" (or "1/1"). Footage is
 * resampled to 1:1 before it is stored, so anything else is refused. */
static int parse_aspect(const char *s) {
    long n, d;
    if (parse_fps(s, &n, &d) != 0) return -1;
    return (n == 1 && d == 1) ? 0 : -1;
}

static uint64_t specs_frame_bytes(const specs_t *sp) {
    const sampling_info_t *si = &SAMPLINGS[sp->sampling];
    uint64_t y = (uint64_t)sp->width * (uint64_t)sp->height;
    uint64_t cw = ((uint64_t)sp->width + (uint64_t)si->hsub - 1) / (uint64_t)si->hsub;
    uint64_t ch = ((uint64_t)sp->height + (uint64_t)si->vsub - 1) / (uint64_t)si->vsub;
    return (y + 2 * cw * ch) * (sp->bitdepth > 8 ? 2u : 1u);
}

/* The y4m stream header line for these specs (with trailing newline).
 * Pixel aspect is always 1:1 (A1:1). Interlacing isn't in specs.txt, so
 * it is written as progressive (Ip). Returns snprintf's result. */
static int y4m_header(const specs_t *sp, char *buf, size_t n) {
    const sampling_info_t *si = &SAMPLINGS[sp->sampling];
    return snprintf(buf, n, "YUV4MPEG2 W%d H%d F%ld:%ld Ip A1:1 C%s\n",
                    sp->width, sp->height, sp->fps_num, sp->fps_den,
                    sp->bitdepth > 8 ? si->y4m10 : si->y4m8);
}

/* Returns 0 = loaded, 1 = file does not exist, -1 = unreadable/invalid
 * (reason in err). An existing but invalid file is never overwritten. */
static int specs_load(const char *path, specs_t *sp, char *err, size_t errlen) {
    FILE *f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT) return 1;
        snprintf(err, errlen, "cannot open: %s", strerror(errno));
        return -1;
    }
    int have_res = 0, have_smp = 0, have_fps = 0, have_bd = 0;
    sp->have_aspect = 0;
    char line[256];
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        char *t = trim(line);
        if (*t == 0 || *t == '#') continue;
        char *eq = strchr(t, '=');
        if (!eq) {
            snprintf(err, errlen, "line %d: expected key=value", lineno);
            fclose(f);
            return -1;
        }
        *eq = 0;
        char *key = trim(t);
        char *val = trim(eq + 1);
        for (char *c = key; *c; c++) *c = (char)tolower((unsigned char)*c);

        if (!strcmp(key, "resolution")) {
            if (have_res) { snprintf(err, errlen, "line %d: duplicate resolution", lineno); fclose(f); return -1; }
            if (parse_resolution(val, &sp->width, &sp->height) != 0) {
                snprintf(err, errlen, "line %d: bad resolution '%s' (expected WxH, e.g. resolution=640x480)", lineno, val);
                fclose(f);
                return -1;
            }
            have_res = 1;
        } else if (!strcmp(key, "sampling")) {
            if (have_smp) { snprintf(err, errlen, "line %d: duplicate sampling", lineno); fclose(f); return -1; }
            sp->sampling = parse_sampling(val);
            if (sp->sampling < 0) {
                snprintf(err, errlen, "line %d: unknown sampling '%s' (valid: 420 422 444)", lineno, val);
                fclose(f);
                return -1;
            }
            have_smp = 1;
        } else if (!strcmp(key, "framerate")) {
            if (have_fps) { snprintf(err, errlen, "line %d: duplicate framerate", lineno); fclose(f); return -1; }
            if (parse_fps(val, &sp->fps_num, &sp->fps_den) != 0) {
                snprintf(err, errlen, "line %d: bad framerate '%s' (must be a fraction N:D, e.g. framerate=24000:1001)", lineno, val);
                fclose(f);
                return -1;
            }
            have_fps = 1;
        } else if (!strcmp(key, "bitdepth")) {
            if (have_bd) { snprintf(err, errlen, "line %d: duplicate bitdepth", lineno); fclose(f); return -1; }
            if (parse_bitdepth(val, &sp->bitdepth) != 0) {
                snprintf(err, errlen, "line %d: bad bitdepth '%s' (valid: 8 10)", lineno, val);
                fclose(f);
                return -1;
            }
            have_bd = 1;
        } else if (!strcmp(key, "aspect")) {
            if (sp->have_aspect) { snprintf(err, errlen, "line %d: duplicate aspect", lineno); fclose(f); return -1; }
            if (parse_aspect(val) != 0) {
                snprintf(err, errlen, "line %d: unsupported aspect '%s' (only aspect=1:1 -- resample to square pixels first)", lineno, val);
                fclose(f);
                return -1;
            }
            sp->have_aspect = 1;
        } else {
            fprintf(stderr, "y4mstore: %s: line %d: ignoring unknown key '%s'\n", path, lineno, key);
        }
    }
    fclose(f);
    if (!have_res || !have_smp || !have_fps || !have_bd) {
        snprintf(err, errlen, "missing %s",
                 !have_res ? "resolution" : !have_smp ? "sampling" :
                 !have_fps ? "framerate" : "bitdepth");
        return -1;
    }
    sp->frame_bytes = specs_frame_bytes(sp);
    return 0;
}

/* Prints "  <key>> ", reads one line, and returns the value (trimmed, with an
 * optional leading "key=" removed), or NULL on EOF. */
static char *prompt_value(const char *key, char *buf, size_t n) {
    fprintf(stderr, "  %s> ", key);
    fflush(stderr);
    if (!fgets(buf, (int)n, stdin)) return NULL;
    return strip_key_prefix(trim(buf), key);
}

/* Ask on the terminal until each answer parses; every question shows an
 * example line and the y4m header token it produces. Returns 0, or -1 if
 * input ends (EOF) before all four values are in. Prompts go to stderr so
 * they stay visible when stdout is redirected or piped. */
static int specs_prompt(specs_t *sp) {
    char buf[128];
    char *v;

    fprintf(stderr, "\n  Resolution: width x height in pixels (y4m W and H).\n"
                    "    example line:  resolution=640x480       -> W640 H480\n");
    for (;;) {
        if (!(v = prompt_value("resolution", buf, sizeof(buf)))) return -1;
        if (parse_resolution(v, &sp->width, &sp->height) == 0) break;
        fprintf(stderr, "    not valid -- expected e.g. resolution=640x480 (each side 1..65535)\n");
    }

    fprintf(stderr, "\n  Sampling: chroma subsampling (y4m C tag). Valid values:\n");
    for (size_t i = 0; i < N_SAMPLINGS; i++) {
        fprintf(stderr, "    %-4s %-40s -> C%s\n", SAMPLINGS[i].name, SAMPLINGS[i].desc, SAMPLINGS[i].y4m8);
    }
    fprintf(stderr, "    frame files must hold planes in y4m order: Y, then U, then V.\n"
                    "    example line:  sampling=420\n");
    for (;;) {
        if (!(v = prompt_value("sampling", buf, sizeof(buf)))) return -1;
        sp->sampling = parse_sampling(v);
        if (sp->sampling >= 0) break;
        fprintf(stderr, "    not valid -- pick one of 420, 422, 444, e.g. sampling=420\n");
    }

    fprintf(stderr, "\n  Frame rate: a fraction N:D, frames per second (y4m F).\n"
                    "    example line:  framerate=24000:1001     -> F24000:1001\n"
                    "    also valid:    24:1   25:1   30000:1001   60:1   60000:1001\n");
    for (;;) {
        if (!(v = prompt_value("framerate", buf, sizeof(buf)))) return -1;
        if (parse_fps(v, &sp->fps_num, &sp->fps_den) == 0) break;
        fprintf(stderr, "    not valid -- must be a fraction N:D, e.g. framerate=24000:1001 "
                        "(write 24 as 24:1)\n");
    }

    fprintf(stderr, "\n  Bit depth: bits per sample. Valid values: 8 10\n"
                    "    8  = one byte per sample                  -> C%s / C%s / C%s\n"
                    "    10 = two bytes per sample, little-endian  -> C%s / C%s / C%s\n"
                    "    example line:  bitdepth=10\n",
            SAMPLINGS[0].y4m8, SAMPLINGS[1].y4m8, SAMPLINGS[2].y4m8,
            SAMPLINGS[0].y4m10, SAMPLINGS[1].y4m10, SAMPLINGS[2].y4m10);
    for (;;) {
        if (!(v = prompt_value("bitdepth", buf, sizeof(buf)))) return -1;
        if (parse_bitdepth(v, &sp->bitdepth) == 0) break;
        fprintf(stderr, "    not valid -- pick 8 or 10, e.g. bitdepth=10\n");
    }

    sp->frame_bytes = specs_frame_bytes(sp);
    return 0;
}

/* Exclusive create ("wx") so an existing file is never clobbered. */
static int specs_write(const char *path, const specs_t *sp) {
    FILE *f = fopen(path, "wx");
    if (!f) return -1;
    fprintf(f, "resolution=%dx%d\nsampling=%s\nframerate=%ld:%ld\nbitdepth=%d\naspect=1:1\n",
            sp->width, sp->height, SAMPLINGS[sp->sampling].name,
            sp->fps_num, sp->fps_den, sp->bitdepth);
    if (fclose(f) != 0) return -1;
    return 0;
}

/* Appends "aspect=1:1" to an existing specs.txt (adding a newline first if
 * the file doesn't end with one). Returns 0 on success. */
static int specs_add_aspect(const char *path) {
    FILE *f = fopen(path, "a+");
    if (!f) return -1;
    int need_nl = 0;
    if (fseek(f, -1, SEEK_END) == 0) need_nl = (fgetc(f) != '\n');
    fseek(f, 0, SEEK_END);   /* required between a read and a write */
    if (fprintf(f, "%saspect=1:1\n", need_nl ? "\n" : "") < 0) { fclose(f); return -1; }
    return fclose(f) == 0 ? 0 : -1;
}

static void specs_print(const specs_t *sp, const char *path) {
    fprintf(stderr, "y4mstore: specs (%s): %dx%d %s %d-bit, %ld:%ld fps, %llu bytes/frame\n",
            path, sp->width, sp->height, SAMPLINGS[sp->sampling].name, sp->bitdepth,
            sp->fps_num, sp->fps_den, (unsigned long long)sp->frame_bytes);
}

/* Find <dir>/specs.txt; if absent, ask the user and write it. Returns 0
 * with sp filled in, or -1 (message already printed). */
static int specs_for_dir(const char *dir, specs_t *sp, int quiet) {
    char path[PATH_MAX];
    snprintf(path, sizeof(path), "%s/%s", dir, SPECS_FILENAME);
    char err[256];
    int rc = specs_load(path, sp, err, sizeof(err));
    if (rc < 0) {
        fprintf(stderr, "y4mstore: %s: %s\n", path, err);
        return -1;
    }
    if (rc == 0) {
        if (!sp->have_aspect) {
            /* Shown even with -q: this edits the user's file. */
            if (specs_add_aspect(path) == 0) {
                fprintf(stderr, "y4mstore: added aspect=1:1 to %s (square pixels; resample to 1:1 before storing)\n", path);
            } else {
                fprintf(stderr, "y4mstore: warning: could not add aspect=1:1 to %s: %s (using 1:1 for this run)\n",
                        path, strerror(errno));
            }
        }
        if (!quiet) specs_print(sp, path);
        return 0;
    }

    if (!isatty(STDIN_FILENO)) {
        fprintf(stderr,
            "y4mstore: %s not found, and stdin is not a terminal so it cannot be asked for.\n"
            "         Create it with these lines (aspect=1:1 is added automatically if omitted):\n"
            "           resolution=640x480\n"
            "           sampling=420\n"
            "           framerate=24000:1001\n"
            "           bitdepth=8\n"
            "           aspect=1:1\n", path);
        return -1;
    }
    fprintf(stderr, "y4mstore: %s not found -- describe the raw YUV in this set\n"
                    "         (applies to every clip subdirectory). Type either the value\n"
                    "         or the whole example line.\n", path);
    if (specs_prompt(sp) != 0) {
        fprintf(stderr, "\ny4mstore: input ended before the specs were complete; nothing written\n");
        return -1;
    }
    if (specs_write(path, sp) != 0) {
        fprintf(stderr, "y4mstore: warning: could not write %s: %s (using these values for this run only)\n",
                path, strerror(errno));
    } else if (!quiet) {
        fprintf(stderr, "y4mstore: wrote %s (aspect=1:1 added automatically)\n", path);
    }
    if (!quiet) specs_print(sp, path);
    return 0;
}

/* Counts the regular files sitting directly in dir, other than specs.txt
 * and dotfiles (symlinks followed). A directory with any of these holds
 * frames, so it is a CLIP directory; a set directory holds only
 * specs.txt and clip subdirectories. */
static size_t count_stray_top_level_files(const char *dir) {
    DIR *d = opendir(dir);
    if (!d) return 0;
    size_t n = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] == '.') continue;
        if (!strcmp(ent->d_name, SPECS_FILENAME)) continue;
        char full[PATH_MAX];
        snprintf(full, sizeof(full), "%s/%s", dir, ent->d_name);
        struct stat st;
        if (stat(full, &st) == 0 && S_ISREG(st.st_mode)) n++;
    }
    closedir(d);
    return n;
}

/* ---------------------------------------------------------------- */
/* Path collection                                                  */
/* ---------------------------------------------------------------- */

typedef struct {
    char **paths;
    size_t count;
    size_t cap;
} pathlist_t;

static void pathlist_init(pathlist_t *pl) {
    pl->paths = NULL;
    pl->count = 0;
    pl->cap = 0;
}

static void pathlist_push(pathlist_t *pl, const char *p) {
    if (pl->count == pl->cap) {
        pl->cap = pl->cap ? pl->cap * 2 : 1024;
        pl->paths = realloc(pl->paths, pl->cap * sizeof(char *));
    }
    pl->paths[pl->count++] = strdup(p);
}

static void scan_dir(const char *dir, pathlist_t *pl, int quiet, int top) {
    DIR *d = opendir(dir);
    if (!d) {
        if (!quiet) fprintf(stderr, "y4mstore: warning: cannot open %s: %s\n",
                             dir, strerror(errno));
        return;
    }
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, "..")) continue;
        char full[PATH_MAX];
        snprintf(full, sizeof(full), "%s/%s", dir, ent->d_name);
        struct stat st;
        if (lstat(full, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) {
            scan_dir(full, pl, quiet, 0);
        } else if (S_ISREG(st.st_mode)) {
            /* specs.txt describes the set; it is not a frame. */
            if (top && !strcmp(ent->d_name, SPECS_FILENAME)) continue;
            pathlist_push(pl, full);
        }
    }
    closedir(d);
}

static void read_stdin_list(pathlist_t *pl) {
    char line[PATH_MAX];
    while (fgets(line, sizeof(line), stdin)) {
        size_t l = strlen(line);
        while (l && (line[l - 1] == '\n' || line[l - 1] == '\r')) line[--l] = 0;
        if (l) pathlist_push(pl, line);
    }
}

/* ---------------------------------------------------------------- */
/* Serve mode (-s): stream a clip directory as yuv4mpeg (y4m)        */
/*                                                                   */
/* The named directory holds the raw frame files (planar Y,U,V at    */
/* the specs' resolution/sampling/bitdepth). specs.txt is read from   */
/* the directory's PARENT (the set directory). Each file becomes one  */
/* or more y4m frames, in natural filename order: header, then        */
/* "FRAME\n" + raw bytes per frame. Raw frames are already in y4m     */
/* plane order, so the bytes are copied through untouched.            */
/* ---------------------------------------------------------------- */

typedef struct {
    char *name;
    uint64_t size;
} frame_file_t;

/* Filename order where digit runs compare as numbers: f2 < f10. */
static int natural_cmp(const char *a, const char *b) {
    while (*a && *b) {
        if (isdigit((unsigned char)*a) && isdigit((unsigned char)*b)) {
            const char *ra = a, *rb = b;              /* whole digit runs */
            while (isdigit((unsigned char)*ra)) ra++;
            while (isdigit((unsigned char)*rb)) rb++;
            const char *za = a, *zb = b;              /* skip leading zeros */
            while (za < ra - 1 && *za == '0') za++;
            while (zb < rb - 1 && *zb == '0') zb++;
            size_t la = (size_t)(ra - za), lb = (size_t)(rb - zb);
            if (la != lb) return la < lb ? -1 : 1;
            int c = memcmp(za, zb, la);
            if (c) return c < 0 ? -1 : 1;
            a = ra;
            b = rb;
        } else {
            if (*a != *b) return (unsigned char)*a < (unsigned char)*b ? -1 : 1;
            a++;
            b++;
        }
    }
    if (*a) return 1;
    if (*b) return -1;
    return 0;
}

static int frame_file_cmp(const void *a, const void *b) {
    const frame_file_t *fa = a, *fb = b;
    int c = natural_cmp(fa->name, fb->name);
    return c ? c : strcmp(fa->name, fb->name);
}

/* Textual parent of a directory path ("set/clip" -> "set", "clip" -> ".",
 * "." -> "./.."). */
static void parent_of(const char *dir, char *out, size_t n) {
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof(tmp), "%s", dir);
    size_t len = strlen(tmp);
    while (len > 1 && tmp[len - 1] == '/') tmp[--len] = 0;
    char *slash = strrchr(tmp, '/');
    const char *base = slash ? slash + 1 : tmp;
    if (!strcmp(base, ".") || !strcmp(base, "..")) {
        snprintf(out, n, "%s/..", tmp);
    } else if (!slash) {
        snprintf(out, n, ".");
    } else if (slash == tmp) {
        snprintf(out, n, "/");
    } else {
        *slash = 0;
        snprintf(out, n, "%s", tmp);
    }
}

static int path_exists(const char *path) {
    struct stat st;
    return stat(path, &st) == 0;
}

/* If path doesn't exist, returns 1 (nothing to confirm) with no output.
 * If it exists and stdin is a terminal, asks the user and returns 1 for
 * yes, 0 (with a message) for no. If it exists and stdin is NOT a
 * terminal, refuses outright (with a message) rather than guessing or
 * hanging on a read that will never come -- same policy as the specs.txt
 * prompt elsewhere in this file. Never appends: the caller always opens
 * the file fresh ("w") after this returns 1, exactly once, and (for -d,
 * which can be written to across several plane sorts in one run) keeps
 * that single handle open rather than reopening it. */
static int confirm_overwrite(const char *path, int quiet) {
    if (!path_exists(path)) return 1;
    if (!isatty(STDIN_FILENO)) {
        fprintf(stderr, "y4mstore: %s already exists, and stdin is not a terminal so it cannot be asked "
                        "whether to overwrite it. Remove it, choose another path, or run interactively.\n", path);
        return 0;
    }
    fprintf(stderr, "y4mstore: %s already exists. Overwrite it? [y/N] ", path);
    fflush(stderr);
    char line[16];
    if (!fgets(line, sizeof(line), stdin) || (line[0] != 'y' && line[0] != 'Y')) {
        if (!quiet) fprintf(stderr, "y4mstore: not overwriting %s; nothing written.\n", path);
        return 0;
    }
    return 1;
}

/* Is there a regular file directly in dir whose name ends in ext (".yuv",
 * ".y", ...; case-insensitive)? */
static int dir_has_ext_files(const char *dir, const char *ext) {
    DIR *d = opendir(dir);
    if (!d) return 0;
    size_t el = strlen(ext);
    struct dirent *ent;
    int found = 0;
    while (!found && (ent = readdir(d)) != NULL) {
        size_t l = strlen(ent->d_name);
        if (ent->d_name[0] != '.' && l > el && strcasecmp(ent->d_name + l - el, ext) == 0) found = 1;
    }
    closedir(d);
    return found;
}

static int has_yuv_files(const char *dir) { return dir_has_ext_files(dir, ".yuv"); }

/* A "3-component" (elementary) directory, the kind -e writes: y/, u/ and
 * v/ subdirectories holding .y, .u and .v files, and no whole-frame .yuv
 * files of its own (a directory with those is a clip, whatever else it
 * contains). The file extensions are what make this specific: a set that
 * merely has clips named y, u and v is not mistaken for one. */
static int is_component_dir(const char *dir) {
    static const char *const PL[3] = {"y", "u", "v"};
    static const char *const EX[3] = {".y", ".u", ".v"};
    char sub[PATH_MAX + 8];
    if (strlen(dir) >= PATH_MAX - 4) return 0;
    for (int p = 0; p < 3; p++) {
        snprintf(sub, sizeof(sub), "%s/%s", dir, PL[p]);
        if (!dir_has_ext_files(sub, EX[p])) return 0;
    }
    return !has_yuv_files(dir);
}

/* True if dir's immediate contents are ONLY subdirectories (plus,
 * optionally, specs.txt) -- no loose regular files directly inside it.
 * That shape means dir is a set directory holding more than one clip
 * (or exactly one, nested), not a clip itself: pointing a sort at it
 * would otherwise recurse blindly through every clip's files as one
 * undifferentiated pile, mixing separate clips (and, for component
 * clips, mixing y/u/v planes together) with no regard for the
 * structure. Callers should refuse rather than silently flattening it. */
static int dir_holds_only_subdirs(const char *dir) {
    DIR *d = opendir(dir);
    if (!d) return 0;
    int any_subdir = 0, any_loose_file = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, "..")) continue;
        char full[PATH_MAX];
        snprintf(full, sizeof(full), "%s/%s", dir, ent->d_name);
        struct stat st;
        if (lstat(full, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) any_subdir = 1;
        else if (S_ISREG(st.st_mode) && strcmp(ent->d_name, SPECS_FILENAME) != 0) any_loose_file = 1;
    }
    closedir(d);
    return any_subdir && !any_loose_file;
}

/* Byte sizes of a frame's y, u and v planes at these specs. */
static void plane_sizes(const specs_t *sp, uint64_t psz[3]) {
    const sampling_info_t *si = &SAMPLINGS[sp->sampling];
    const uint64_t bps = sp->bitdepth > 8 ? 2 : 1;
    psz[0] = (uint64_t)sp->width * (uint64_t)sp->height * bps;
    psz[1] = psz[2] = (((uint64_t)sp->width + (uint64_t)si->hsub - 1) / (uint64_t)si->hsub) *
                      (((uint64_t)sp->height + (uint64_t)si->vsub - 1) / (uint64_t)si->vsub) * bps;
}

static int name_natural_cmp(const void *a, const void *b) {
    const char *x = *(char *const *)a, *y = *(char *const *)b;
    int c = natural_cmp(x, y);
    return c ? c : strcmp(x, y);
}

/* The regular files in dir whose names end in ext (or, with ext == NULL,
 * every regular file except dotfiles and specs.txt: a clip's frame files), as "dir/name" strings in
 * natural filename order (so the run does not depend on directory order).
 * Returns a malloc'd array (or NULL when empty); *count is its length. */
static char **list_ext_files(const char *dir, const char *ext, size_t *count) {
    DIR *d = opendir(dir);
    *count = 0;
    if (!d) return NULL;
    size_t el = ext ? strlen(ext) : 0, cap = 0, n = 0;
    char **names = NULL;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        size_t l = strlen(ent->d_name);
        if (ent->d_name[0] == '.') continue;
        if (ext) {
            if (l <= el || strcasecmp(ent->d_name + l - el, ext) != 0) continue;
        } else if (!strcmp(ent->d_name, SPECS_FILENAME)) {
            continue;
        }
        char full[PATH_MAX + 300];
        snprintf(full, sizeof(full), "%.*s/%s", PATH_MAX - 1, dir, ent->d_name);
        struct stat st;
        if (stat(full, &st) != 0 || !S_ISREG(st.st_mode)) continue;
        if (n == cap) {
            cap = cap ? cap * 2 : 1024;
            char **nn = realloc(names, cap * sizeof(*names));
            if (!nn) { closedir(d); for (size_t i = 0; i < n; i++) free(names[i]); free(names); return NULL; }
            names = nn;
        }
        names[n] = strdup(full);
        if (!names[n]) { closedir(d); for (size_t i = 0; i < n; i++) free(names[i]); free(names); return NULL; }
        n++;
    }
    closedir(d);
    if (n) qsort(names, n, sizeof(*names), name_natural_cmp);
    *count = n;
    return names;
}


/* Which directory's specs.txt applies to the directory the user named?
 *
 *  - A directory holding no frame files is a SET directory: its own
 *    specs.txt applies.
 *  - A directory holding frame files, or y/u/v component planes (made by
 *    -e), is a CLIP-like directory: the specs.txt
 *    belongs to its parent (the set). A copy inside the clip is ignored
 *    with a warning; if the parent has none, the clip's own copy is used
 *    (with a note) rather than asking again; if neither exists the
 *    parent is where a new one will be created.
 *
 * Decided by what the directory contains, never by whether a specs.txt
 * happens to sit inside it. Writes the directory to use into out. */
static void specs_locate(const char *dir, char *out, size_t n, int quiet) {
    if (count_stray_top_level_files(dir) == 0 && !is_component_dir(dir)) {
        snprintf(out, n, "%s", dir);
        return;
    }
    char own[PATH_MAX + 32], up[PATH_MAX + 32], parent[PATH_MAX + 8];
    parent_of(dir, parent, sizeof(parent));
    {   /* "." / ".." read badly in messages ("move it up to ."): spell those out */
        size_t pl = strlen(parent);
        if (!strcmp(parent, ".") || !strcmp(parent, "..") ||
            (pl >= 3 && !strcmp(parent + pl - 3, "/.."))) {
            char abs[PATH_MAX];
            if (realpath(parent, abs)) snprintf(parent, sizeof(parent), "%s", abs);
        }
    }
    snprintf(own, sizeof(own), "%s/%s", dir, SPECS_FILENAME);
    snprintf(up, sizeof(up), "%s/%s", parent, SPECS_FILENAME);
    int own_has = path_exists(own);
    if (path_exists(up)) {
        if (own_has && !quiet) {
            fprintf(stderr, "y4mstore: warning: ignoring %s -- a clip directory's specs come from its set "
                            "directory (%s); this copy can be deleted\n", own, up);
        }
        snprintf(out, n, "%s", parent);
    } else if (own_has) {
        if (!quiet) {
            fprintf(stderr, "y4mstore: note: no %s; using %s instead. Move it up to %s to match the layout.\n",
                    up, own, parent);
        }
        snprintf(out, n, "%s", dir);
    } else {
        snprintf(out, n, "%s", parent);
    }
}

#define SERVE_BUF (1u << 20)

/* Returns the process exit status. A closed pipe (player quit) is a
 * normal stop, not an error. */
static int serve_write_failed(void) {
    if (errno == EPIPE) return 0;
    fprintf(stderr, "\ny4mstore: write error: %s\n", strerror(errno));
    return 1;
}

/* Text quoted from an input stream (or a filename) is untrusted: a damaged or
 * hostile source must not be able to put control characters or terminal
 * escape sequences in a message. Copies at most n-1 bytes, showing anything that
 * isn't printable ASCII as '?'. */
static const char *safe_text(const char *in, char *out, size_t n) {
    size_t i = 0;
    for (; in[i] && i + 1 < n; i++) {
        unsigned char c = (unsigned char)in[i];
        out[i] = (c >= 0x20 && c < 0x7f) ? (char)c : '?';
    }
    out[i] = 0;
    return out;
}

/* ---- serving a 3-component directory ------------------------------------ */
/* The reverse of -e: a directory of y/, u/ and v/ planes goes back out as a */
/* y4m stream, each frame being its y, u and v files one after another.     */
/* The three planes of a frame are paired by FILENAME STEM (yuv-0000007.y,  */
/* .u and .v), never by position, so a missing or extra file can't quietly  */
/* mix planes of different frames; if any frame lacks a plane, or a plane   */
/* is the wrong size, nothing is written. Frames go out in natural stem     */
/* order, exactly as whole-frame -s orders its files.                       */

typedef struct {
    char *stem;
    char *path;
} pfile_t;

static int stem_cmp(const char *a, const char *b) {
    int c = natural_cmp(a, b);
    return c ? c : strcmp(a, b);
}

static int pfile_cmp(const void *a, const void *b) {
    return stem_cmp(((const pfile_t *)a)->stem, ((const pfile_t *)b)->stem);
}

/* "dir/yuv-0000001.y" with extlen 2 -> "yuv-0000001" (malloc'd) */
static char *plane_stem(const char *path, size_t extlen) {
    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    size_t l = strlen(base);
    if (l < extlen) return NULL;
    char *st = malloc(l - extlen + 1);
    if (!st) return NULL;
    memcpy(st, base, l - extlen);
    st[l - extlen] = 0;
    return st;
}

static int serve_planes(const char *dir, const char *out_path, int quiet) {
    static const char PLANE[3] = {'y', 'u', 'v'};
    static const char *const EXT[3] = {".y", ".u", ".v"};
    int rc = 1;
    pfile_t *pf[3] = {NULL, NULL, NULL};
    size_t np[3] = {0, 0, 0};
    size_t (*fr)[3] = NULL;          /* complete frames: an index into each pf[p] */
    size_t nfr = 0;
    FILE *out = NULL;
    unsigned char *buf = NULL;

    /* 1. the three plane lists, each in stem order */
    for (int p = 0; p < 3; p++) {
        char sub[PATH_MAX + 8];
        snprintf(sub, sizeof(sub), "%.*s/%c", PATH_MAX - 4, dir, PLANE[p]);
        char **paths = list_ext_files(sub, EXT[p], &np[p]);
        if (!paths || np[p] == 0) { fprintf(stderr, "y4mstore: no %s files found in %s\n", EXT[p], sub); goto done; }
        pf[p] = calloc(np[p], sizeof(pfile_t));
        if (!pf[p]) { for (size_t i = 0; i < np[p]; i++) free(paths[i]); free(paths); fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
        for (size_t i = 0; i < np[p]; i++) {
            pf[p][i].path = paths[i];
            pf[p][i].stem = plane_stem(paths[i], strlen(EXT[p]));
            if (!pf[p][i].stem) { fprintf(stderr, "y4mstore: out of memory\n"); free(paths); goto done; }
        }
        free(paths);
        qsort(pf[p], np[p], sizeof(pfile_t), pfile_cmp);
    }

    /* 2. the specs: a component directory's come from its parent, like a clip's */
    char sdir[PATH_MAX + 8];
    specs_locate(dir, sdir, sizeof(sdir), quiet);
    specs_t sp;
    if (specs_for_dir(sdir, &sp, quiet) != 0) goto done;
    uint64_t psz[3];
    plane_sizes(&sp, psz);

    /* 3. pair the planes by stem. Every frame must have all three. */
    size_t maxfr = np[0] > np[1] ? np[0] : np[1];
    if (np[2] > maxfr) maxfr = np[2];
    fr = malloc(maxfr * sizeof(*fr));
    if (!fr) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
    {
        size_t i[3] = {0, 0, 0}, incomplete = 0;
        while (i[0] < np[0] || i[1] < np[1] || i[2] < np[2]) {
            const char *minst = NULL;
            for (int p = 0; p < 3; p++) {
                if (i[p] < np[p] && (!minst || stem_cmp(pf[p][i[p]].stem, minst) < 0)) minst = pf[p][i[p]].stem;
            }
            int has[3], nhas = 0;
            for (int p = 0; p < 3; p++) {
                has[p] = i[p] < np[p] && stem_cmp(pf[p][i[p]].stem, minst) == 0;
                nhas += has[p];
            }
            if (nhas == 3) {
                fr[nfr][0] = i[0]; fr[nfr][1] = i[1]; fr[nfr][2] = i[2];
                nfr++;
            } else {
                if (incomplete < 5) {
                    char found[8] = "", missing[8] = "";
                    for (int p = 0; p < 3; p++) {
                        size_t l = strlen(has[p] ? found : missing);
                        char *dst = has[p] ? found : missing;
                        if (l) dst[l++] = ' ';
                        dst[l++] = PLANE[p];
                        dst[l] = 0;
                    }
                    char shown[96];
                    fprintf(stderr, "y4mstore: frame '%s': found %s, missing %s\n", safe_text(minst, shown, sizeof(shown)), found, missing);
                }
                incomplete++;
            }
            for (int p = 0; p < 3; p++) if (has[p]) i[p]++;
        }
        if (incomplete) {
            fprintf(stderr, "y4mstore: %zu frame name%s not present in all of y/, u/ and v/ (y=%zu, u=%zu, v=%zu files). "
                            "Nothing was written.\n", incomplete, incomplete == 1 ? " is" : "s are", np[0], np[1], np[2]);
            goto done;
        }
    }

    /* 4. every plane file must be exactly the plane size these specs give */
    {
        size_t bad = 0;
        for (size_t f = 0; f < nfr; f++) {
            for (int p = 0; p < 3; p++) {
                struct stat st;
                const char *path = pf[p][fr[f][p]].path;
                if (stat(path, &st) != 0 || (uint64_t)st.st_size != psz[p]) {
                    if (bad < 5) {
                        if (stat(path, &st) != 0) fprintf(stderr, "y4mstore: %s: %s\n", path, strerror(errno));
                        else fprintf(stderr, "y4mstore: %s: %llu bytes, expected %llu (one %c plane)\n", path,
                                     (unsigned long long)st.st_size, (unsigned long long)psz[p], PLANE[p]);
                    }
                    bad++;
                }
            }
        }
        if (bad) {
            fprintf(stderr, "y4mstore: %zu plane file%s do not match the specs (%dx%d %s %d-bit: planes of %llu, %llu, %llu bytes). "
                            "Nothing was written. Check %s/%s.\n", bad, bad == 1 ? "" : "s", sp.width, sp.height,
                    SAMPLINGS[sp.sampling].name, sp.bitdepth, (unsigned long long)psz[0], (unsigned long long)psz[1],
                    (unsigned long long)psz[2], sdir, SPECS_FILENAME);
            goto done;
        }
    }

    /* 5. the stream */
    char hdr[256];
    int hl = y4m_header(&sp, hdr, sizeof(hdr));
    if (hl < 0 || (size_t)hl >= sizeof(hdr)) { fprintf(stderr, "y4mstore: internal error building header\n"); goto done; }
    if (!quiet) {
        fprintf(stderr, "y4mstore: stream header: %.*s\n", hl - 1, hdr);
        fprintf(stderr, "y4mstore: serving %zu frame%s from y/u/v planes (%zu files, %.1f s at %ld:%ld fps)\n",
                nfr, nfr == 1 ? "" : "s", 3 * nfr, (double)nfr * (double)sp.fps_den / (double)sp.fps_num,
                sp.fps_num, sp.fps_den);
    }
    out = stdout;
    if (out_path) {
        out = fopen(out_path, "wb");
        if (!out) { fprintf(stderr, "y4mstore: cannot open %s for writing: %s\n", out_path, strerror(errno)); out = NULL; goto done; }
    }
    signal(SIGPIPE, SIG_IGN);
    setvbuf(out, NULL, _IOFBF, SERVE_BUF);
    buf = malloc(SERVE_BUF);
    if (!buf) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
    if (fwrite(hdr, 1, (size_t)hl, out) != (size_t)hl) { rc = serve_write_failed(); goto done; }

    for (size_t f = 0; f < nfr; f++) {
        if (fwrite("FRAME\n", 1, 6, out) != 6) { rc = serve_write_failed(); goto done; }
        for (int p = 0; p < 3; p++) {
            const char *path = pf[p][fr[f][p]].path;
            FILE *in = fopen(path, "rb");
            if (!in) { fprintf(stderr, "\ny4mstore: cannot read %s: %s\n", path, strerror(errno)); goto done; }
            uint64_t left = psz[p];
            while (left) {
                size_t want = left < SERVE_BUF ? (size_t)left : SERVE_BUF;
                size_t r = fread(buf, 1, want, in);
                if (r == 0) { fprintf(stderr, "\ny4mstore: %s ended early (file changed while serving?)\n", path); fclose(in); goto done; }
                if (fwrite(buf, 1, r, out) != r) { fclose(in); rc = serve_write_failed(); goto done; }
                left -= r;
            }
            fclose(in);
        }
        if (out_path && !quiet) progress_bar("serving", f + 1, nfr);
    }
    if (fflush(out) != 0) { rc = serve_write_failed(); goto done; }
    rc = 0;
    if (out_path && !quiet) fprintf(stderr, "y4mstore: wrote %zu frame%s to %s\n", nfr, nfr == 1 ? "" : "s", out_path);

done:
    if (out && out != stdout) {
        if (fclose(out) != 0 && rc == 0) { fprintf(stderr, "y4mstore: error closing %s: %s\n", out_path, strerror(errno)); rc = 1; }
    }
    free(buf);
    free(fr);
    for (int p = 0; p < 3; p++) {
        if (pf[p]) {
            for (size_t i = 0; i < np[p]; i++) { free(pf[p][i].stem); free(pf[p][i].path); }
            free(pf[p]);
        }
    }
    return rc;
}

static int serve_clip(const char *dir, const char *out_path, int quiet) {
    if (!out_path && isatty(STDOUT_FILENO)) {
        fprintf(stderr, "y4mstore: refusing to write a video stream to a terminal.\n"
                        "         Pipe it (y4mstore -s %s | ffplay -), redirect it (> clip.y4m), or use -f FILE.\n", dir);
        return 1;
    }

    /* y/, u/ and v/ planes (what -e writes) go back out as a stream too */
    if (is_component_dir(dir)) return serve_planes(dir, out_path, quiet);

    /* Frame files: regular files directly in dir, skipping dotfiles and
     * any specs.txt (which describes frames, it isn't one). */
    DIR *d = opendir(dir);
    if (!d) {
        fprintf(stderr, "y4mstore: cannot open %s: %s\n", dir, strerror(errno));
        return 1;
    }
    frame_file_t *files = NULL;
    size_t nfiles = 0, cap = 0, subdirs = 0;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] == '.') continue;
        if (!strcmp(ent->d_name, SPECS_FILENAME)) continue;
        char full[PATH_MAX];
        snprintf(full, sizeof(full), "%s/%s", dir, ent->d_name);
        struct stat st;
        if (stat(full, &st) != 0) continue;
        if (S_ISDIR(st.st_mode)) { subdirs++; continue; }
        if (!S_ISREG(st.st_mode)) continue;
        if (nfiles == cap) {
            cap = cap ? cap * 2 : 1024;
            files = realloc(files, cap * sizeof(*files));
            if (!files) { closedir(d); fprintf(stderr, "y4mstore: out of memory\n"); return 1; }
        }
        files[nfiles].name = strdup(ent->d_name);
        files[nfiles].size = (uint64_t)st.st_size;
        nfiles++;
    }
    closedir(d);

    int rc = 1;
    FILE *out = NULL;
    unsigned char *buf = NULL;

    if (nfiles == 0) {
        char own[PATH_MAX + 32];
        snprintf(own, sizeof(own), "%s/%s", dir, SPECS_FILENAME);
        if (subdirs || path_exists(own)) {
            fprintf(stderr, "y4mstore: %s has no frame files of its own, so it looks like a set directory.\n"
                            "         Name one of its clip subdirectories, e.g. y4mstore -s %s/<clip>\n", dir, dir);
        } else {
            fprintf(stderr, "y4mstore: no frame files in %s\n", dir);
        }
        goto done;
    }
    if (subdirs && !quiet) {
        fprintf(stderr, "y4mstore: warning: ignoring %zu subdirector%s in %s "
                        "(frames are read from the top level only)\n",
                subdirs, subdirs == 1 ? "y" : "ies", dir);
    }
    qsort(files, nfiles, sizeof(*files), frame_file_cmp);

    char sdir[PATH_MAX + 8];   /* the directory whose specs.txt applies */
    specs_locate(dir, sdir, sizeof(sdir), quiet);
    specs_t sp;
    if (specs_for_dir(sdir, &sp, quiet) != 0) goto done;

    /* Check every file before writing a byte: each must hold a whole
     * number of frames at these specs. */
    size_t bad = 0;
    uint64_t frames = 0;
    for (size_t i = 0; i < nfiles; i++) {
        if (files[i].size == 0 || files[i].size % sp.frame_bytes != 0) {
            if (bad < 5) {
                fprintf(stderr, "y4mstore: %s/%s: %llu bytes is not a multiple of the %llu-byte frame size\n",
                        dir, files[i].name, (unsigned long long)files[i].size,
                        (unsigned long long)sp.frame_bytes);
            }
            bad++;
        } else {
            frames += files[i].size / sp.frame_bytes;
        }
    }
    if (bad) {
        fprintf(stderr, "y4mstore: %zu of %zu file%s do not match the specs (%dx%d %s %d-bit = %llu bytes/frame).\n"
                        "         Nothing was written. Check %s/%s.\n",
                bad, nfiles, nfiles == 1 ? "" : "s", sp.width, sp.height,
                SAMPLINGS[sp.sampling].name, sp.bitdepth,
                (unsigned long long)sp.frame_bytes, sdir, SPECS_FILENAME);
        goto done;
    }

    char hdr[256];
    int hl = y4m_header(&sp, hdr, sizeof(hdr));
    if (hl < 0 || (size_t)hl >= sizeof(hdr)) { fprintf(stderr, "y4mstore: internal error building header\n"); goto done; }

    if (!quiet) {
        fprintf(stderr, "y4mstore: stream header: %.*s\n", hl - 1, hdr);
        fprintf(stderr, "y4mstore: serving %llu frame%s from %zu file%s (%.1f s at %ld:%ld fps)\n",
                (unsigned long long)frames, frames == 1 ? "" : "s", nfiles, nfiles == 1 ? "" : "s",
                (double)frames * (double)sp.fps_den / (double)sp.fps_num, sp.fps_num, sp.fps_den);
    }

    out = stdout;
    if (out_path) {
        out = fopen(out_path, "wb");
        if (!out) {
            fprintf(stderr, "y4mstore: cannot open %s for writing: %s\n", out_path, strerror(errno));
            out = NULL;
            goto done;
        }
    }
    signal(SIGPIPE, SIG_IGN);
    setvbuf(out, NULL, _IOFBF, SERVE_BUF);
    buf = malloc(SERVE_BUF);
    if (!buf) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }

    if (fwrite(hdr, 1, (size_t)hl, out) != (size_t)hl) { rc = serve_write_failed(); goto done; }

    uint64_t sent = 0;
    for (size_t i = 0; i < nfiles; i++) {
        char full[PATH_MAX];
        snprintf(full, sizeof(full), "%s/%s", dir, files[i].name);
        FILE *in = fopen(full, "rb");
        if (!in) {
            fprintf(stderr, "\ny4mstore: cannot read %s: %s\n", full, strerror(errno));
            goto done;
        }
        uint64_t nf = files[i].size / sp.frame_bytes;
        for (uint64_t k = 0; k < nf; k++) {
            if (fwrite("FRAME\n", 1, 6, out) != 6) { fclose(in); rc = serve_write_failed(); goto done; }
            uint64_t left = sp.frame_bytes;
            while (left) {
                size_t want = left < SERVE_BUF ? (size_t)left : SERVE_BUF;
                size_t r = fread(buf, 1, want, in);
                if (r == 0) {
                    fprintf(stderr, "\ny4mstore: %s ended early (file changed while serving?)\n", full);
                    fclose(in);
                    goto done;
                }
                if (fwrite(buf, 1, r, out) != r) { fclose(in); rc = serve_write_failed(); goto done; }
                left -= r;
            }
            sent++;
            if (out_path && !quiet) progress_bar("serving", (size_t)sent, (size_t)frames);
        }
        fclose(in);
    }
    if (fflush(out) != 0) { rc = serve_write_failed(); goto done; }
    rc = 0;
    if (out_path && !quiet) {
        fprintf(stderr, "y4mstore: wrote %llu frame%s to %s\n",
                (unsigned long long)frames, frames == 1 ? "" : "s", out_path);
    }

done:
    if (out && out != stdout) {
        if (fclose(out) != 0 && rc == 0) {
            fprintf(stderr, "y4mstore: error closing %s: %s\n", out_path, strerror(errno));
            rc = 1;
        }
    }
    free(buf);
    for (size_t i = 0; i < nfiles; i++) free(files[i].name);
    free(files);
    return rc;
}

/* ---------------------------------------------------------------- */
/* Alternate hash (-a): additive component averages                  */
/*                                                                   */
/* No similarity hash and no O(n^2) passes. Each frame file (exactly */
/* one frame at the specs' resolution/sampling/bitdepth) is split    */
/* into its Y, U and V planes. Every plane's samples are added up    */
/* and divided by that plane's own sample count (345600 for the Y    */
/* plane of a 720x480 frame; 4:2:0 chroma planes have 86400).        */
/* Each average is multiplied by 100 before rounding (ADDITIVE_SCALE),  */
/* since an average of many samples has real, usable precision far     */
/* finer than one sample's own bit depth -- see the ADDITIVE_SCALE     */
/* comment above alt_key_digits() for why 100 specifically. This makes */
/* each plane's digit width bit-depth-dependent (5 digits at 8-bit,    */
/* up to 25500; 6 digits at 10-bit, up to 102300), and the three are   */
/* joined into one number, Y first:                                    */
/*     8-bit, Y=250.00 U=127.00 V=36.00  ->  250001270003600           */
/* Frames are then sorted numerically by that number. One sort.        */
/*                                                                   */
/* Equal keys do NOT mean identical frames: moving a pixel changes   */
/* nothing about an average. So ties are ordered by filename (which  */
/* keeps runs deterministic) and are never treated as duplicates.    */
/* ---------------------------------------------------------------- */

typedef struct {
    uint64_t key;        /* Y/U/V packed as one decimal number, digit width per alt_key_digits() */
    const char *path;
} alt_entry_t;

/* Precision multiplier for -a's averaged key. An average of many samples
 * has real, usable resolution far finer than one sample's own bit depth
 * (its noise floor shrinks as 1/sqrt(sample count)); *1 (the original
 * scheme) rounded to a single sample's own integer range and threw that
 * away, especially on the U/V planes (4x fewer samples than Y, so 4x
 * coarser floor, using the SAME 4-digit range as Y). *100 sits solidly
 * past a conservative *10 step and comfortably below where real, non-CG
 * source noise would start to dominate -- see conversation for the
 * quantization/noise-floor analysis this was chosen from. */
#define ADDITIVE_SCALE 100

/* Decimal digits needed to hold (max sample value at this bit depth) *
 * ADDITIVE_SCALE -- e.g. 255*100=25500 (5 digits) at 8-bit, 1023*100=
 * 102300 (6 digits) at 10-bit. This replaces the old fixed width of 4,
 * which only happened to fit both bit depths at *1. */
static int alt_key_digits(int bitdepth) {
    uint64_t v = (uint64_t)(bitdepth > 8 ? 1023 : 255) * ADDITIVE_SCALE;
    int d = 1;
    while (v >= 10) { v /= 10; d++; }
    return d;
}

/* 10^n as a uint64_t, n small (at most ~6 here). */
static uint64_t upow10(int n) {
    uint64_t v = 1;
    while (n-- > 0) v *= 10;
    return v;
}

/* Rounded, scaled average of n samples; 8-bit samples are bytes, wider
 * ones are 16-bit little-endian words. Clamped to max_val (the true
 * maximum possible average at this bit depth and scale) as a safety net;
 * a correct average can never exceed it. */
static uint64_t alt_plane_avg(const unsigned char *p, uint64_t n, int wide, uint64_t max_val) {
    uint64_t sum = 0;
    if (!wide) {
        for (uint64_t i = 0; i < n; i++) sum += p[i];
    } else {
        for (uint64_t i = 0; i < n; i++) sum += (uint64_t)p[2 * i] | ((uint64_t)p[2 * i + 1] << 8);
    }
    uint64_t avg = (sum * ADDITIVE_SCALE + n / 2) / n;
    return avg > max_val ? max_val : avg;
}

static uint64_t alt_frame_key(const unsigned char *frame, const specs_t *sp) {
    const sampling_info_t *si = &SAMPLINGS[sp->sampling];
    int wide = sp->bitdepth > 8;
    uint64_t bps = wide ? 2 : 1;
    uint64_t ys = (uint64_t)sp->width * (uint64_t)sp->height;
    uint64_t cw = ((uint64_t)sp->width + (uint64_t)si->hsub - 1) / (uint64_t)si->hsub;
    uint64_t ch = ((uint64_t)sp->height + (uint64_t)si->vsub - 1) / (uint64_t)si->vsub;
    uint64_t cs = cw * ch;
    int d = alt_key_digits(sp->bitdepth);
    uint64_t max_val = (uint64_t)(wide ? 1023 : 255) * ADDITIVE_SCALE;
    uint64_t y = alt_plane_avg(frame, ys, wide, max_val);
    uint64_t u = alt_plane_avg(frame + ys * bps, cs, wide, max_val);
    uint64_t v = alt_plane_avg(frame + (ys + cs) * bps, cs, wide, max_val);
    return y * upow10(2 * d) + u * upow10(d) + v;
}

typedef struct {
    char **paths;
    size_t n;
    const specs_t *sp;
    uint64_t *keys;
    int *err;                 /* 0 = ok, otherwise an errno value */
    atomic_size_t *next;      /* next frame to claim */
    atomic_size_t *done;
    uint64_t expect;          /* exact size of each file, in bytes */
    int plane_only;           /* 1: each file is ONE plane; its key is that plane's average */
} alt_job_t;

/* Per-frame parallelism: every frame is independent, so each worker
 * claims the next single frame the moment it finishes its last one --
 * balanced load with no thread-per-frame spawn cost. */
static void *alt_worker(void *arg) {
    alt_job_t *j = (alt_job_t *)arg;
    unsigned char *buf = malloc((size_t)j->expect);
    for (;;) {
        size_t i = atomic_fetch_add(j->next, 1);
        if (i >= j->n) break;
        if (!buf) {
            j->err[i] = ENOMEM;
        } else {
            FILE *f = fopen(j->paths[i], "rb");
            if (!f) {
                j->err[i] = errno ? errno : EIO;
            } else {
                if (fread(buf, 1, (size_t)j->expect, f) != (size_t)j->expect) {
                    j->err[i] = EIO;
                } else if (j->plane_only) {
                    int wide = j->sp->bitdepth > 8;
                    j->keys[i] = alt_plane_avg(buf, j->expect / (wide ? 2 : 1), wide,
                                               (uint64_t)(wide ? 1023 : 255) * ADDITIVE_SCALE);
                } else {
                    j->keys[i] = alt_frame_key(buf, j->sp);
                }
                fclose(f);
            }
        }
        atomic_fetch_add(j->done, 1);
    }
    free(buf);
    return NULL;
}

static int alt_entry_cmp(const void *a, const void *b) {
    const alt_entry_t *ea = a, *eb = b;
    if (ea->key != eb->key) return ea->key < eb->key ? -1 : 1;
    int c = natural_cmp(ea->path, eb->path);
    return c ? c : strcmp(ea->path, eb->path);
}

/* Returns the process exit status. Nothing is written unless every file
 * is exactly the expected size and readable.
 *
 *  - plane_bytes == 0: each file is one whole frame (the -a mode); the key
 *    is the Y/U/V number (digit width per alt_key_digits(), x3).
 *  - plane_bytes  > 0: each file is ONE plane of exactly that many bytes
 *    (used by -e on the y, u and v directories); the key is that plane's
 *    single average, digit width per alt_key_digits().
 *  - out_shared != NULL: write the sorted paths there (the caller owns the
 *    stream); otherwise open out_path, or use stdout.
 *  - label, if given, prefixes the summary lines ("y", "u", "v"). */
static int alt_sort(char **paths, size_t n, const specs_t *sp, const char *specs_dir,
                    const char *out_path, FILE *out_shared, uint64_t plane_bytes, const char *label,
                    int nthreads, int quiet, FILE *dump) {
    const uint64_t expect = plane_bytes ? plane_bytes : sp->frame_bytes;
    const int digits = plane_bytes ? alt_key_digits(sp->bitdepth) : 3 * alt_key_digits(sp->bitdepth);
    char lab[16] = "";
    if (label) snprintf(lab, sizeof(lab), "[%s] ", label);
    int rc = 1;
    uint64_t *keys = NULL;
    int *err = NULL;
    alt_entry_t *ent = NULL;
    FILE *out = NULL;

    /* Every file must be exactly one frame at these specs. */
    size_t bad = 0;
    for (size_t i = 0; i < n; i++) {
        struct stat st;
        if (stat(paths[i], &st) != 0 || (uint64_t)st.st_size != expect) {
            if (bad < 5) {
                if (stat(paths[i], &st) != 0)
                    fprintf(stderr, "y4mstore: %s: %s\n", paths[i], strerror(errno));
                else
                    fprintf(stderr, "y4mstore: %s: %llu bytes, expected %llu (one %s)\n", paths[i],
                            (unsigned long long)st.st_size, (unsigned long long)expect, plane_bytes ? "plane" : "frame");
            }
            bad++;
        }
    }
    if (bad) {
        fprintf(stderr, "y4mstore: %zu of %zu file%s are not exactly one %s (%dx%d %s %d-bit = %llu bytes).\n"
                        "         -a needs one %s per file. Nothing was written. Check %s/%s.\n",
                bad, n, n == 1 ? "" : "s", plane_bytes ? "plane" : "frame", sp->width, sp->height,
                SAMPLINGS[sp->sampling].name, sp->bitdepth, (unsigned long long)expect,
                plane_bytes ? "plane" : "frame", specs_dir, SPECS_FILENAME);
        return 1;
    }

    keys = calloc(n, sizeof(*keys));
    err = calloc(n, sizeof(*err));
    ent = malloc(n * sizeof(*ent));
    if (!keys || !err || !ent) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }

    int nt = nthreads < 1 ? 1 : nthreads;
    if ((size_t)nt > n) nt = (int)n;
    atomic_size_t next = 0, finished = 0;
    alt_job_t job = { paths, n, sp, keys, err, &next, &finished, expect, plane_bytes != 0 };
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)nt);
    if (!th) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    int launched = 0;
    for (int t = 0; t < nt; t++) {
        if (pthread_create(&th[t], NULL, alt_worker, &job) == 0) launched++;
    }
    if (launched == 0) { fprintf(stderr, "y4mstore: cannot start worker threads\n"); free(th); goto done; }
    if (!quiet) {
        size_t seen;
        do {
            struct timespec nap = {0, 50 * 1000 * 1000};
            nanosleep(&nap, NULL);
            seen = atomic_load(&finished);
            progress_bar("averaging", seen, n);
        } while (seen < n);
    }
    for (int t = 0; t < launched; t++) pthread_join(th[t], NULL);
    free(th);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    size_t failed = 0;
    for (size_t i = 0; i < n; i++) {
        if (err[i]) {
            if (failed < 5) fprintf(stderr, "y4mstore: cannot read %s: %s\n", paths[i], strerror(err[i]));
            failed++;
        }
    }
    if (failed) {
        fprintf(stderr, "y4mstore: %zu file%s could not be read. Nothing was written.\n", failed, failed == 1 ? "" : "s");
        goto done;
    }

    for (size_t i = 0; i < n; i++) { ent[i].key = keys[i]; ent[i].path = paths[i]; }
    qsort(ent, n, sizeof(*ent), alt_entry_cmp);

    if (!quiet) {
        size_t distinct = 0, shared = 0;
        for (size_t i = 0; i < n; ) {
            size_t j = i + 1;
            while (j < n && ent[j].key == ent[i].key) j++;
            distinct++;
            if (j - i > 1) shared += j - i;
            i = j;
        }
        double secs = (double)(t1.tv_sec - t0.tv_sec) + (double)(t1.tv_nsec - t0.tv_nsec) / 1e9;
        const char *unit = plane_bytes ? "plane" : "frame";
        fprintf(stderr, "y4mstore: %saveraged %zu %s%s (%llu bytes each) in %.2f s with %d thread%s\n",
                lab, n, unit, n == 1 ? "" : "s", (unsigned long long)expect, secs, launched, launched == 1 ? "" : "s");
        fprintf(stderr, "y4mstore: %s%zu distinct key%s, lowest %0*llu, highest %0*llu; %zu %s%s share a key "
                        "(ordered by filename; equal keys are not duplicates)\n",
                lab, distinct, distinct == 1 ? "" : "s", digits, (unsigned long long)ent[0].key,
                digits, (unsigned long long)ent[n - 1].key, shared, unit, shared == 1 ? "" : "s");
    }
    if (dump) {
        for (size_t i = 0; i < n; i++) fprintf(dump, "%0*llu  %s\n", digits, (unsigned long long)ent[i].key, ent[i].path);
        fflush(dump);
    }

    if (out_shared) {
        out = out_shared;               /* the caller's stream: written to, never closed here */
    } else {
        out = stdout;
        if (out_path) {
            out = fopen(out_path, "w");
            if (!out) {
                fprintf(stderr, "y4mstore: cannot open %s for writing: %s\n", out_path, strerror(errno));
                out = NULL;
                goto done;
            }
        }
    }
    for (size_t i = 0; i < n; i++) fprintf(out, "%s\n", ent[i].path);
    fflush(out);
    rc = 0;
    if (!out_shared && out_path && !quiet) fprintf(stderr, "y4mstore: wrote %zu file path(s) to %s\n", n, out_path);

done:
    if (out && out != stdout && out != out_shared) {
        if (fclose(out) != 0 && rc == 0) {
            fprintf(stderr, "y4mstore: error closing %s: %s\n", out_path, strerror(errno));
            rc = 1;
        }
    }
    free(keys);
    free(err);
    free(ent);
    return rc;
}

/* ---------------------------------------------------------------- */
/* Ingest (-i): a y4m stream becomes a clip directory                */
/*                                                                   */
/*   y4mstore -i FILE|- <clipdir>                                     */
/*                                                                   */
/* The stream header is turned into specs (resolution, sampling,     */
/* framerate, bitdepth). specs.txt in the clip's parent (the set) is */
/* created from it if absent; if present it must be identical, and   */
/* otherwise nothing at all is written. Each frame then becomes one  */
/* raw file, yiv-000001.yuv onward, in the new clip directory. The   */
/* frame bytes are copied through untouched -- y4m frames already    */
/* are raw planar Y,U,V in exactly the layout the other modes read.  */
/* ---------------------------------------------------------------- */

#define INGEST_LINE_MAX 4096
#define INGEST_BUF (1u << 20)

/* y4m colorspace tags this tool can store, and the specs they map to. */
static const struct {
    const char *tag, *sampling;
    int depth;
    const char *note;
} Y4M_COLORSPACES[] = {
    {"420jpeg",  "420",  8, NULL},
    {"420",      "420",  8, NULL},
    {"420mpeg2", "420",  8, "420mpeg2 chroma siting is not stored; the clip is recorded as plain 420"},
    {"420paldv", "420",  8, "420paldv chroma siting is not stored; the clip is recorded as plain 420"},
    {"422",      "422",  8, NULL},
    {"444",      "444",  8, NULL},
    {"420p10",   "420", 10, NULL},
    {"422p10",   "422", 10, NULL},
    {"444p10",   "444", 10, NULL},
};

static void notes_add(char *notes, size_t n, const char *fmt, ...) __attribute__((format(printf, 3, 4)));
static void notes_add(char *notes, size_t n, const char *fmt, ...) {
    size_t l = strlen(notes);
    if (l >= n) return;
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(notes + l, n - l, fmt, ap);
    va_end(ap);
}

/* Reads one line (without its newline). Returns 0 = got a line,
 * 1 = clean EOF before any byte, -1 = EOF in the middle of a line,
 * -2 = line longer than the buffer. */
static int ingest_read_line(FILE *in, char *line, size_t n) {
    size_t l = 0;
    int c = getc(in);
    if (c == EOF) return 1;
    for (;;) {
        if (c == '\n') { line[l] = 0; return 0; }
        if (c == EOF) return -1;
        if (l + 1 >= n) return -2;
        line[l++] = (char)c;
        c = getc(in);
    }
}

static int parse_int_strict(const char *v, long lo, long hi, long *out) {
    char *e;
    long x = strtol(v, &e, 10);
    if (e == v || *e != 0 || x < lo || x > hi) return -1;
    *out = x;
    return 0;
}

/* Turn a y4m stream header line into specs. Returns 0, or -1 with the
 * reason in err. Anything the tool can't store faithfully is refused
 * rather than guessed at. Informational remarks go into notes. */
static int y4m_parse_header(const char *line, specs_t *sp, char *err, size_t errlen,
                            char *notes, size_t noteslen) {
    char buf[INGEST_LINE_MAX];
    char shown[80];   /* sanitized copy of a stream token, for messages */
    snprintf(buf, sizeof(buf), "%s", line);
    int have_w = 0, have_h = 0, have_f = 0;
    const char *cs = "420jpeg";   /* the y4m default when there is no C tag */
    char interlace = '?';
    char *save = NULL;
    char *tok = strtok_r(buf, " ", &save);
    if (!tok || strcmp(tok, "YUV4MPEG2") != 0) {
        snprintf(err, errlen, "not a yuv4mpeg stream (it does not start with YUV4MPEG2)");
        return -1;
    }
    while ((tok = strtok_r(NULL, " ", &save)) != NULL) {
        const char *v = tok + 1;
        long x;
        switch (tok[0]) {
        case 'W':
            if (parse_int_strict(v, 1, 65535, &x) != 0) { snprintf(err, errlen, "bad width tag '%s'", safe_text(tok, shown, sizeof(shown))); return -1; }
            sp->width = (int)x; have_w = 1;
            break;
        case 'H':
            if (parse_int_strict(v, 1, 65535, &x) != 0) { snprintf(err, errlen, "bad height tag '%s'", safe_text(tok, shown, sizeof(shown))); return -1; }
            sp->height = (int)x; have_h = 1;
            break;
        case 'F':
            if (parse_fps(v, &sp->fps_num, &sp->fps_den) != 0) {
                snprintf(err, errlen, "the stream has no usable frame rate ('%s'); specs.txt needs a real one",
                         safe_text(tok, shown, sizeof(shown)));
                return -1;
            }
            have_f = 1;
            break;
        case 'I':
            interlace = v[0] ? v[0] : '?';
            break;
        case 'A': {
            char *e;
            long n = strtol(v, &e, 10);
            long d = 0;
            if (e == v || *e != ':' || (d = strtol(e + 1, &e, 10), *e != 0) || n < 0 || d < 0) {
                snprintf(err, errlen, "bad pixel aspect tag '%s'", safe_text(tok, shown, sizeof(shown)));
                return -1;
            }
            if (n == 0 && d == 0) {
                notes_add(notes, noteslen, "y4mstore: note: the stream's pixel aspect is unknown (A0:0); treated as square (1:1)\n");
            } else if (n == 0 || d == 0 || n != d) {
                snprintf(err, errlen, "the stream has non-square pixels (%s); resample to square pixels (1:1) before ingesting",
                         safe_text(tok, shown, sizeof(shown)));
                return -1;
            }
            break;
        }
        case 'C':
            cs = v;
            break;
        default:
            break;   /* X comments and tags this tool doesn't use */
        }
    }
    if (!have_w || !have_h) { snprintf(err, errlen, "the header has no %s", !have_w ? "width (W)" : "height (H)"); return -1; }
    if (!have_f) { snprintf(err, errlen, "the header has no frame rate (F tag); specs.txt needs one"); return -1; }

    int found = -1;
    for (size_t i = 0; i < sizeof(Y4M_COLORSPACES) / sizeof(Y4M_COLORSPACES[0]); i++) {
        if (!strcmp(cs, Y4M_COLORSPACES[i].tag)) { found = (int)i; break; }
    }
    if (found < 0) {
        snprintf(err, errlen, "unsupported colorspace 'C%s'. Supported: 420 (420jpeg/420mpeg2/420paldv), 422, 444, "
                              "and their 10-bit forms 420p10, 422p10, 444p10", safe_text(cs, shown, sizeof(shown)));
        return -1;
    }
    sp->sampling = parse_sampling(Y4M_COLORSPACES[found].sampling);
    sp->bitdepth = Y4M_COLORSPACES[found].depth;
    sp->have_aspect = 1;
    sp->frame_bytes = specs_frame_bytes(sp);
    if (Y4M_COLORSPACES[found].note) notes_add(notes, noteslen, "y4mstore: note: %s\n", Y4M_COLORSPACES[found].note);
    if (interlace == 't' || interlace == 'b' || interlace == 'm') {
        notes_add(notes, noteslen, "y4mstore: note: the stream is interlaced (I%c); frames are stored as they come "
                                   "and specs.txt does not record interlacing\n", interlace);
    }
    return 0;
}

/* Creates <dir>/specs.txt with its complete content in a single atomic
 * step, and fails with EEXIST if it already exists. (Plain exclusive
 * create makes the file empty first and fills it a moment later, so a
 * second process running at the same time could read a half-written
 * file.) The temp file is a dotfile so it can never look like a frame or
 * a stray file if the process is killed mid-way. Filesystems without hard
 * links fall back to the plain exclusive create. */
static int specs_write_atomic(const char *dir, const specs_t *sp) {
    char path[PATH_MAX + 32], tmp[PATH_MAX + 64];
    snprintf(path, sizeof(path), "%s/%s", dir, SPECS_FILENAME);
    snprintf(tmp, sizeof(tmp), "%s/.%s.tmp%ld", dir, SPECS_FILENAME, (long)getpid());
    unlink(tmp);   /* a stale one from a killed run with the same pid */
    if (specs_write(tmp, sp) != 0) return -1;
    int r = link(tmp, path);
    int e = errno;
    unlink(tmp);
    if (r == 0) return 0;
    if (e == EEXIST) { errno = EEXIST; return -1; }
    return specs_write(path, sp);
}

/* Print each spec that differs between specs.txt (a) and the stream (b). */
static size_t specs_diff(const specs_t *a, const specs_t *b) {
    size_t n = 0;
    char x[64], y[64];
    if (a->width != b->width || a->height != b->height) {
        snprintf(x, sizeof(x), "%dx%d", a->width, a->height);
        snprintf(y, sizeof(y), "%dx%d", b->width, b->height);
        fprintf(stderr, "  %-11s specs.txt %-14s stream %s\n", "resolution", x, y);
        n++;
    }
    if (a->sampling != b->sampling) {
        fprintf(stderr, "  %-11s specs.txt %-14s stream %s\n", "sampling", SAMPLINGS[a->sampling].name, SAMPLINGS[b->sampling].name);
        n++;
    }
    if (a->bitdepth != b->bitdepth) {
        snprintf(x, sizeof(x), "%d", a->bitdepth);
        snprintf(y, sizeof(y), "%d", b->bitdepth);
        fprintf(stderr, "  %-11s specs.txt %-14s stream %s\n", "bitdepth", x, y);
        n++;
    }
    if (a->fps_num != b->fps_num || a->fps_den != b->fps_den) {
        snprintf(x, sizeof(x), "%ld:%ld", a->fps_num, a->fps_den);
        snprintf(y, sizeof(y), "%ld:%ld", b->fps_num, b->fps_den);
        fprintf(stderr, "  %-11s specs.txt %-14s stream %s\n", "framerate", x, y);
        n++;
    }
    return n;
}

/* ---- live frame counter (-p) --------------------------------------------- */
/* A piped stream has no known size, so there is no progress bar; without   */
/* this a long "ffmpeg ... | y4mstore -i -" prints nothing until it is done.  */
/* -p shows a running count of frames read. On a terminal it is one line    */
/* redrawn about 4 times a second; if stderr is redirected (a log file) it  */
/* is a plain line every couple of seconds instead, so a log doesn't fill   */
/* with carriage-return redraws. It is opt-in because upstream tools such   */
/* as ffmpeg draw their own \r status line on the same terminal, and two    */
/* programs redrawing one line garble each other. Only the reading thread   */
/* ever touches it.                                                         */

static int g_progress_count = 0;      /* -p was given */

static struct {
    int on;          /* counting for this run */
    int tty;         /* stderr is a terminal */
    int drawn;       /* a \r line is on screen: end it before printing anything else */
    double t0, last;
} g_ctr;

static double mono_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

static void counter_start(int enable) {
    memset(&g_ctr, 0, sizeof(g_ctr));
    g_ctr.on = enable;
    g_ctr.tty = isatty(STDERR_FILENO);
    g_ctr.t0 = g_ctr.last = mono_now();
}

static void counter_draw(const char *label, uint64_t frames, uint64_t bytes) {
    double now = mono_now(), el = now - g_ctr.t0;
    double fps = el > 0.05 ? (double)frames / el : 0.0;
    int secs = (int)el;
    fprintf(stderr, g_ctr.tty ? "\ry4mstore: %s: %llu frames, %.1f MB, %.1f frames/s, %d:%02d elapsed   "
                              : "y4mstore: %s: %llu frames, %.1f MB, %.1f frames/s, %d:%02d elapsed\n",
            label, (unsigned long long)frames, (double)bytes / 1e6, fps, secs / 60, secs % 60);
    g_ctr.drawn = g_ctr.tty;
    g_ctr.last = now;
}

/* Called after every frame; draws only when enough time has passed. */
static void counter_update(const char *label, uint64_t frames, uint64_t bytes) {
    if (!g_ctr.on) return;
    if (mono_now() - g_ctr.last < (g_ctr.tty ? 0.25 : 2.0)) return;
    counter_draw(label, frames, bytes);
}

/* Ends a redrawn line, so the next message starts on a fresh one. */
static void counter_break(void) {
    if (g_ctr.drawn) {
        fputc('\n', stderr);
        g_ctr.drawn = 0;
    }
}

/* The last state, then a clean line end. (In a log the summary line that
 * follows already says the total, so only a terminal gets a final redraw.) */
static void counter_finish(const char *label, uint64_t frames, uint64_t bytes) {
    if (!g_ctr.on) return;
    if (g_ctr.tty) counter_draw(label, frames, bytes);
    counter_break();
}

/* Opens the y4m source for -i / -e: "-" is stdin (refused if it is a
 * terminal), anything else is a file (a directory is refused with the
 * way to feed one in). flag is "-i" or "-e", for the messages. */
static int y4m_open_source(const char *src, const char *flag, const char *dest, FILE **in_out, off_t *fsize_out) {
    FILE *in;
    if (!strcmp(src, "-")) {
        if (isatty(STDIN_FILENO)) {
            fprintf(stderr, "y4mstore: refusing to read a video stream from a terminal; pipe it in "
                            "(ffmpeg ... -f yuv4mpegpipe - | y4mstore %s - %s)\n", flag, dest);
            return 1;
        }
        in = stdin;
    } else {
        struct stat dst;
        if (stat(src, &dst) == 0 && S_ISDIR(dst.st_mode)) {
            fprintf(stderr, "y4mstore: %s is a directory; %s reads a y4m file or stream.\n"
                            "         To start from a clip directory, pipe another instance in: "
                            "y4mstore -s %s | y4mstore %s - %s\n", src, flag, src, flag, dest);
            return 1;
        }
        in = fopen(src, "rb");
        if (!in) {
            fprintf(stderr, "y4mstore: cannot open %s: %s\n", src, strerror(errno));
            return 1;
        }
    }
    off_t fsize = 0;
    {
        struct stat st;
        if (in != stdin && fstat(fileno(in), &st) == 0 && S_ISREG(st.st_mode)) fsize = st.st_size;
    }
    setvbuf(in, NULL, _IOFBF, INGEST_BUF);
    *in_out = in;
    *fsize_out = fsize;
    return 0;
}

/* Reads and parses the stream header into hs, and says what the stream
 * is. hdr_bytes gets the header's length including its newline. */
static int y4m_read_header(FILE *in, specs_t *hs, int quiet, uint64_t *hdr_bytes) {
    char line[INGEST_LINE_MAX], err[512], notes[1024] = "";
    int lr = ingest_read_line(in, line, sizeof(line));
    if (lr == 1) { fprintf(stderr, "y4mstore: the input is empty; expected a yuv4mpeg stream\n"); return 1; }
    if (lr != 0) {
        fprintf(stderr, "y4mstore: could not read a yuv4mpeg header line (%s)\n",
                lr == -2 ? "line too long" : "input ended early");
        return 1;
    }
    memset(hs, 0, sizeof(*hs));
    if (y4m_parse_header(line, hs, err, sizeof(err), notes, sizeof(notes)) != 0) {
        fprintf(stderr, "y4mstore: %s\n", err);
        return 1;
    }
    if (!quiet) {
        fprintf(stderr, "y4mstore: stream: %dx%d %s %d-bit, %ld:%ld fps, %llu bytes/frame\n",
                hs->width, hs->height, SAMPLINGS[hs->sampling].name, hs->bitdepth, hs->fps_num, hs->fps_den,
                (unsigned long long)hs->frame_bytes);
        fputs(notes, stderr);
    }
    *hdr_bytes = (uint64_t)strlen(line) + 1;
    return 0;
}

/* specs.txt in the set directory (parent): create it from the stream
 * header, or require that the existing one is identical. Returns 0 to
 * carry on, 1 to stop (message printed, nothing written). */
static int specs_settle(const char *parent, const specs_t *hs, int quiet) {
    char spath[PATH_MAX + 32], err[512];
    snprintf(spath, sizeof(spath), "%s/%s", parent, SPECS_FILENAME);
    for (int attempt = 0; attempt < 2; attempt++) {
        specs_t fs;
        int sr = specs_load(spath, &fs, err, sizeof(err));
        if (sr < 0) {
            fprintf(stderr, "y4mstore: %s: %s\n         The existing specs.txt is left as it is; nothing was written.\n", spath, err);
            return 1;
        }
        if (sr == 0) {
            if (fs.width != hs->width || fs.height != hs->height || fs.sampling != hs->sampling ||
                fs.bitdepth != hs->bitdepth || fs.fps_num != hs->fps_num || fs.fps_den != hs->fps_den) {
                fprintf(stderr, "y4mstore: the stream does not match %s -- stopping, nothing was written:\n", spath);
                specs_diff(&fs, hs);
                fprintf(stderr, "         Use a set whose specs match the stream, or a new set directory.\n");
                return 1;
            }
            if (!fs.have_aspect) {
                if (specs_add_aspect(spath) == 0)
                    fprintf(stderr, "y4mstore: added aspect=1:1 to %s (square pixels; resample to 1:1 before storing)\n", spath);
                else
                    fprintf(stderr, "y4mstore: warning: could not add aspect=1:1 to %s: %s\n", spath, strerror(errno));
            }
            if (!quiet) fprintf(stderr, "y4mstore: %s matches the stream\n", spath);
            return 0;
        }
        if (!path_exists(parent) && mkdir(parent, 0777) != 0 && errno != EEXIST) {
            fprintf(stderr, "y4mstore: cannot create %s: %s\n", parent, strerror(errno));
            return 1;
        }
        if (specs_write_atomic(parent, hs) == 0) {
            if (!quiet) fprintf(stderr, "y4mstore: wrote %s from the stream header (aspect=1:1 added automatically)\n", spath);
            return 0;
        }
        if (errno != EEXIST) {
            fprintf(stderr, "y4mstore: cannot write %s: %s\n", spath, strerror(errno));
            return 1;
        }
        /* EEXIST: another run created it a moment ago -- look again and compare */
    }
    fprintf(stderr, "y4mstore: could not settle %s\n", spath);
    return 1;
}

/* dest must be a new or empty directory. *to_create is set to 1 if it
 * does not exist yet. Returns 0 to carry on, 1 to stop. */
static int dest_new_or_empty(const char *dest, int *to_create, const char *hint) {
    struct stat st;
    *to_create = 0;
    if (stat(dest, &st) == 0) {
        if (!S_ISDIR(st.st_mode)) {
            fprintf(stderr, "y4mstore: %s exists and is not a directory\n", dest);
            return 1;
        }
        DIR *d = opendir(dest);
        if (!d) { fprintf(stderr, "y4mstore: cannot open %s: %s\n", dest, strerror(errno)); return 1; }
        struct dirent *ent;
        size_t entries = 0;
        while ((ent = readdir(d)) != NULL) {
            if (strcmp(ent->d_name, ".") && strcmp(ent->d_name, "..")) entries++;
        }
        closedir(d);
        if (entries) {
            fprintf(stderr, "y4mstore: %s already holds %zu entr%s. %s; nothing was written.\n",
                    dest, entries, entries == 1 ? "y" : "ies", hint);
            return 1;
        }
        return 0;
    }
    if (errno != ENOENT) {
        fprintf(stderr, "y4mstore: cannot stat %s: %s\n", dest, strerror(errno));
        return 1;
    }
    *to_create = 1;
    return 0;
}

static int ingest_y4m(const char *src, const char *dest, int quiet) {
    int rc = 1;
    FILE *in = NULL;
    unsigned char *buf = NULL;
    int created_dest = 0;
    uint64_t frames = 0;
    char last_name[PATH_MAX + 32] = "";
    char line[INGEST_LINE_MAX];

    if (!strcmp(dest, "-")) {
        fprintf(stderr, "y4mstore: -i needs a clip directory to write into\n");
        return 1;
    }
    off_t fsize = 0;
    if (y4m_open_source(src, "-i", dest, &in, &fsize) != 0) return 1;

    /* 1. the header, before anything is touched */
    specs_t hs;
    uint64_t consumed = 0;
    if (y4m_read_header(in, &hs, quiet, &consumed) != 0) goto done;
    if (!quiet) fprintf(stderr, "y4mstore: writing frames to %s\n", dest);

    /* 2. the destination must be new or empty */
    if (dest_new_or_empty(dest, &created_dest, "-i writes only into a new or empty clip directory") != 0) goto done;

    /* 3. specs.txt in the set directory (the clip's parent) */
    char parent[PATH_MAX + 8];
    parent_of(dest, parent, sizeof(parent));
    if (specs_settle(parent, &hs, quiet) != 0) goto done;

    /* 4. the frames */
    if (created_dest && mkdir(dest, 0777) != 0 && errno != EEXIST) {
        fprintf(stderr, "y4mstore: cannot create %s: %s\n", dest, strerror(errno));
        created_dest = 0;
        goto done;
    }
    buf = malloc(INGEST_BUF);
    if (!buf) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }

    counter_start(g_progress_count && !quiet && fsize == 0);   /* a file source has the progress bar */
    if (g_ctr.on) counter_draw("ingesting", 0, 0);              /* visible before the first frame arrives */
    for (;;) {
        int fr = ingest_read_line(in, line, sizeof(line));
        if (fr == 1) break;   /* clean end of stream, on a frame boundary */
        if (fr != 0) {
            counter_break();
            fprintf(stderr, "y4mstore: stream damaged after %llu frame%s: %s\n", (unsigned long long)frames,
                    frames == 1 ? "" : "s", fr == -2 ? "frame header too long" : "input ended inside a frame header");
            goto done;
        }
        if (strncmp(line, "FRAME", 5) != 0 || (line[5] != 0 && line[5] != ' ')) {
            char shown[24];
            counter_break();
            fprintf(stderr, "y4mstore: stream damaged after %llu frame%s: expected FRAME, found '%s'\n",
                    (unsigned long long)frames, frames == 1 ? "" : "s", safe_text(line, shown, sizeof(shown)));
            goto done;
        }
        consumed += (uint64_t)strlen(line) + 1;
        snprintf(last_name, sizeof(last_name), "%s/yiv-%06llu.yuv", dest, (unsigned long long)(frames + 1));
        FILE *o = fopen(last_name, "wb");
        if (!o) { counter_break(); fprintf(stderr, "y4mstore: cannot create %s: %s\n", last_name, strerror(errno)); goto done; }
        uint64_t left = hs.frame_bytes;
        while (left) {
            size_t want = left < INGEST_BUF ? (size_t)left : INGEST_BUF;
            size_t r = fread(buf, 1, want, in);
            if (r == 0) {
                counter_break();
                fprintf(stderr, "y4mstore: stream ended inside frame %llu (%llu of %llu bytes); "
                                "%llu complete frame%s kept in %s\n", (unsigned long long)(frames + 1),
                        (unsigned long long)(hs.frame_bytes - left), (unsigned long long)hs.frame_bytes,
                        (unsigned long long)frames, frames == 1 ? "" : "s", dest);
                fclose(o);
                unlink(last_name);
                goto done;
            }
            if (fwrite(buf, 1, r, o) != r) {
                counter_break();
                fprintf(stderr, "y4mstore: write error on %s: %s\n", last_name, strerror(errno));
                fclose(o);
                unlink(last_name);
                goto done;
            }
            left -= r;
        }
        if (fclose(o) != 0) {
            counter_break();
            fprintf(stderr, "y4mstore: write error on %s: %s\n", last_name, strerror(errno));
            unlink(last_name);
            goto done;
        }
        consumed += hs.frame_bytes;
        frames++;
        if (fsize > 0 && !quiet) progress_bar("ingesting", (size_t)consumed, (size_t)fsize);
        counter_update("ingesting", frames, consumed);
    }
    counter_finish("ingesting", frames, consumed);
    if (fsize > 0 && !quiet) progress_bar("ingesting", (size_t)fsize, (size_t)fsize);
    if (frames == 0) {
        fprintf(stderr, "y4mstore: the stream has a header but no frames\n");
        goto done;
    }
    if (!quiet) {
        fprintf(stderr, "y4mstore: ingested %llu frame%s (%llu bytes each, %.1f s at %ld:%ld fps) into %s\n",
                (unsigned long long)frames, frames == 1 ? "" : "s", (unsigned long long)hs.frame_bytes,
                (double)frames * (double)hs.fps_den / (double)hs.fps_num, hs.fps_num, hs.fps_den, dest);
    }
    rc = 0;

done:
    counter_break();
    /* a directory this run created but never filled is not left behind */
    if (rc != 0 && created_dest && frames == 0) rmdir(dest);
    free(buf);
    if (in && in != stdin) fclose(in);
    return rc;
}

/* ---------------------------------------------------------------- */
/* Check mode (-c): create specs.txt if missing, else list it        */
/* ---------------------------------------------------------------- */

/* "%.3f" with trailing zeros dropped: 24000/1001 -> "23.976", 24/1 -> "24". */
static void fps_text(const specs_t *sp, char *buf, size_t n) {
    snprintf(buf, n, "%.3f", (double)sp->fps_num / (double)sp->fps_den);
    size_t l = strlen(buf);
    while (l > 0 && buf[l - 1] == '0') buf[--l] = 0;
    if (l > 0 && buf[l - 1] == '.') buf[--l] = 0;
}

/* The report -c exists for, so it goes to stdout and is shown even
 * with -q. Each spec is shown with the y4m header token it produces. */
static void specs_list(const specs_t *sp, const char *path) {
    const sampling_info_t *si = &SAMPLINGS[sp->sampling];
    char hdr[256];
    int hl = y4m_header(sp, hdr, sizeof(hdr));
    if (hl < 0 || (size_t)hl >= sizeof(hdr)) hl = 1;
    char fps[32];
    fps_text(sp, fps, sizeof(fps));
    char val[64], tok[64], note[96];

    printf("%s\n", path);

    snprintf(val, sizeof(val), "%dx%d", sp->width, sp->height);
    snprintf(tok, sizeof(tok), "W%d H%d", sp->width, sp->height);
    printf("  %-11s %-14s -> %s\n", "resolution", val, tok);

    snprintf(tok, sizeof(tok), "C%s", sp->bitdepth > 8 ? si->y4m10 : si->y4m8);
    printf("  %-11s %-14s -> %-12s (%s)\n", "sampling", si->name, tok, si->desc);

    snprintf(val, sizeof(val), "%ld:%ld", sp->fps_num, sp->fps_den);
    snprintf(tok, sizeof(tok), "F%ld:%ld", sp->fps_num, sp->fps_den);
    printf("  %-11s %-14s -> %-12s (%s fps)\n", "framerate", val, tok, fps);

    snprintf(val, sizeof(val), "%d", sp->bitdepth);
    if (sp->bitdepth > 8) snprintf(note, sizeof(note), "2 bytes per sample, little-endian (p10 suffix)");
    else                  snprintf(note, sizeof(note), "1 byte per sample (no suffix)");
    printf("  %-11s %-14s -> %s\n", "bitdepth", val, note);

    printf("  %-11s %-14s -> %s\n", "aspect", "1:1", "A1:1");
    printf("  %-11s %llu bytes\n", "frame size", (unsigned long long)sp->frame_bytes);
    printf("  %-11s %.*s\n", "y4m header", hl - 1, hdr);
    fflush(stdout);
}

/* -c: find the specs.txt that applies to the named directory (see
 * specs_locate: a set directory's own, or a clip directory's parent's).
 * Missing -> ask for the values and create it; present -> list it.
 * Nothing is scanned, sorted or served. Exit status 0 = specs.txt exists
 * and is valid. */
static int check_specs(const char *dir, int quiet) {
    char sdir[PATH_MAX + 8];
    specs_locate(dir, sdir, sizeof(sdir), quiet);
    if (strcmp(sdir, dir) != 0 && !quiet) {
        fprintf(stderr, "y4mstore: %s holds %s, so it is a clip directory; "
                        "checking its set directory %s\n", dir,
                is_component_dir(dir) ? "y/u/v component planes" : "frame files", sdir);
    }

    char path[PATH_MAX + 32];
    snprintf(path, sizeof(path), "%s/%s", sdir, SPECS_FILENAME);
    int existed = path_exists(path);

    specs_t sp;
    if (specs_for_dir(sdir, &sp, 1) != 0) return 1;
    if (!existed) {
        if (!path_exists(path)) {
            fprintf(stderr, "y4mstore: %s could not be created\n", path);
            return 1;
        }
        if (!quiet) fprintf(stderr, "y4mstore: created %s (aspect=1:1 added automatically)\n", path);
    }
    specs_list(&sp, path);
    return 0;
}

/* ---------------------------------------------------------------- */
/* The nilsimsa similarity sort of a list of files                   */
/* ---------------------------------------------------------------- */

/* Hashes every file, orders them (greedy nearest-neighbour, then 2-opt),
 * and writes the ordered paths to out. Returns 0, or 1 on failure. The
 * caller owns paths and out. Used by the default mode and, once per
 * plane, by -e. */
static int nilsimsa_sort_paths(char **paths, size_t count, FILE *out,
                               int max_sweeps, int nthreads, int quiet) {
    unsigned char *digests = malloc(count * DIGEST_BYTES);
    if (!digests) {
        fprintf(stderr, "y4mstore: out of memory allocating digest buffer\n");
        return 1;
    }
    uint64_t *dup_hash_all = malloc(count * sizeof(uint64_t));
    uint64_t *dup_len_all = malloc(count * sizeof(uint64_t));
    if (!dup_hash_all || !dup_len_all) {
        fprintf(stderr, "y4mstore: out of memory allocating duplicate-tracking buffer\n");
        free(digests); free(dup_hash_all); free(dup_len_all);
        return 1;
    }

    if (!quiet) fprintf(stderr, "y4mstore: hashing %zu files...\n", count);
    hash_all(paths, count, digests, dup_hash_all, dup_len_all, nthreads, quiet);

    if (!quiet) fprintf(stderr, "y4mstore: ordering...\n");
    size_t *order = greedy_order(digests, count, quiet);
    optimize_order(digests, order, count, max_sweeps, quiet);
    for (size_t i = 0; i < count; i++) {
        fprintf(out, "%s\n", paths[order[i]]);
    }
    fflush(out);
    free(order);

    report_duplicates(dup_hash_all, dup_len_all, count, quiet);

    free(digests);
    free(dup_hash_all);
    free(dup_len_all);
    return 0;
}

/* ---------------------------------------------------------------- */
/* Elementary planes (-e): a y4m stream split into y, u and v files   */
/*                                                                   */
/*   y4mstore -e FILE|- <subdirectory>                                */
/*                                                                   */
/* Every frame of the stream is cut into its three planes:           */
/*     <subdirectory>/y/yuv-0000001.y                                */
/*     <subdirectory>/u/yuv-0000001.u                                */
/*     <subdirectory>/v/yuv-0000001.v                                */
/* If <subdirectory> already exists and holds .yuv files (a whole-   */
/* frame clip, say from -i), the planes go to <subdirectory>-        */
/* elementary instead. The chosen sorter (nilsimsa, or the additive  */
/* averages with -a) then runs on the y, u and v directories in that */
/* order and the three sorted lists are appended into one list.      */
/* ---------------------------------------------------------------- */

/* ---- the threaded split ------------------------------------------------ */
/* The stream can only be read in order, so the calling thread reads one   */
/* whole frame at a time; each frame is then handed to one of the -t       */
/* workers, which writes its y, u and v files. That is per-frame           */
/* parallelism, the same as -a: frames are independent, and a worker takes */
/* the next frame the moment it finishes its last. A frame is only handed  */
/* over once it has been read completely, so a truncated stream can never  */
/* leave a frame with only some of its planes. Memory is bounded: a fixed  */
/* small set of frame buffers is recycled, so the reader waits when the    */
/* workers fall behind.                                                    */

typedef struct {
    unsigned char *data;
    size_t cap;
} split_buf_t;

typedef struct {
    split_buf_t *buf;
    uint64_t index;           /* 1-based frame number */
} split_job_t;

typedef struct {
    pthread_mutex_t mu;
    pthread_cond_t cv_work;   /* a job is waiting, or the reader is done */
    pthread_cond_t cv_free;   /* a buffer came back, or a worker failed */
    split_job_t *queue;
    size_t qcap, qhead, qcount;
    split_buf_t **freelist;
    size_t nfree;
    split_buf_t *bufs;
    size_t nbufs;
    pthread_t *threads;
    int nthreads;
    int shutdown;             /* the reader has finished: workers exit when the queue is empty */
    int failed;               /* a worker could not write a frame: the reader must stop */
    char errmsg[700];
    char (*plane_dir)[PATH_MAX + 64];
    char plane_ch[3];
    uint64_t psz[3];
} split_ctx_t;

/* Writes one frame's three plane files. On failure removes whatever it
 * wrote of this frame and returns -1 with the reason in err. */
static int split_write_frame(const split_ctx_t *c, const unsigned char *data, uint64_t idx, char *err, size_t errlen) {
    char path[3][PATH_MAX + 128];
    uint64_t off = 0;
    for (int p = 0; p < 3; p++) {
        snprintf(path[p], sizeof(path[p]), "%.*s/yuv-%07llu.%c", PATH_MAX + 63, c->plane_dir[p],
                 (unsigned long long)idx, c->plane_ch[p]);
    }
    for (int p = 0; p < 3; p++) {
        FILE *o = fopen(path[p], "wb");
        int ok = (o != NULL);
        int e = errno;
        if (ok && fwrite(data + off, 1, (size_t)c->psz[p], o) != (size_t)c->psz[p]) { ok = 0; e = errno; }
        if (o && fclose(o) != 0 && ok) { ok = 0; e = errno; }
        if (!ok) {
            snprintf(err, errlen, "write error on %.400s: %.120s", path[p], strerror(e));
            for (int q = 0; q <= p; q++) unlink(path[q]);
            return -1;
        }
        off += c->psz[p];
    }
    return 0;
}

static void *split_worker(void *arg) {
    split_ctx_t *c = (split_ctx_t *)arg;
    for (;;) {
        pthread_mutex_lock(&c->mu);
        while (c->qcount == 0 && !c->shutdown) pthread_cond_wait(&c->cv_work, &c->mu);
        if (c->qcount == 0) { pthread_mutex_unlock(&c->mu); break; }   /* done, and nothing left */
        split_job_t job = c->queue[c->qhead];
        c->qhead = (c->qhead + 1) % c->qcap;
        c->qcount--;
        int skip = c->failed;   /* after a failure, stop writing and just hand buffers back */
        pthread_mutex_unlock(&c->mu);

        char err[600] = "";
        int rc = 0;
        if (!skip) rc = split_write_frame(c, job.buf->data, job.index, err, sizeof(err));

        pthread_mutex_lock(&c->mu);
        if (rc != 0 && !c->failed) {
            c->failed = 1;
            snprintf(c->errmsg, sizeof(c->errmsg), "%s (frame %llu)", err, (unsigned long long)job.index);
        }
        c->freelist[c->nfree++] = job.buf;
        pthread_cond_broadcast(&c->cv_free);
        pthread_mutex_unlock(&c->mu);
    }
    return NULL;
}

static int split_pool_start(split_ctx_t *c, int nthreads, size_t nbufs, char (*plane_dir)[PATH_MAX + 64],
                            const char plane_ch[3], const uint64_t psz[3]) {
    memset(c, 0, sizeof(*c));
    pthread_mutex_init(&c->mu, NULL);
    pthread_cond_init(&c->cv_work, NULL);
    pthread_cond_init(&c->cv_free, NULL);
    c->qcap = nbufs;
    c->nbufs = nbufs;
    c->queue = calloc(nbufs, sizeof(*c->queue));
    c->freelist = calloc(nbufs, sizeof(*c->freelist));
    c->bufs = calloc(nbufs, sizeof(*c->bufs));
    c->threads = calloc((size_t)nthreads, sizeof(*c->threads));
    c->plane_dir = plane_dir;
    memcpy(c->plane_ch, plane_ch, 3);
    memcpy(c->psz, psz, sizeof(c->psz));
    if (!c->queue || !c->freelist || !c->bufs || !c->threads) return -1;
    for (size_t i = 0; i < nbufs; i++) c->freelist[c->nfree++] = &c->bufs[i];
    for (int t = 0; t < nthreads; t++) {
        if (pthread_create(&c->threads[c->nthreads], NULL, split_worker, c) == 0) c->nthreads++;
    }
    return c->nthreads > 0 ? 0 : -1;
}

/* Tells the workers no more frames are coming, waits for them to finish
 * everything already handed over, and frees the pool. */
static void split_pool_finish(split_ctx_t *c) {
    pthread_mutex_lock(&c->mu);
    c->shutdown = 1;
    pthread_cond_broadcast(&c->cv_work);
    pthread_mutex_unlock(&c->mu);
    for (int t = 0; t < c->nthreads; t++) pthread_join(c->threads[t], NULL);
    for (size_t i = 0; i < c->nbufs; i++) free(c->bufs[i].data);
    free(c->queue);
    free(c->freelist);
    free(c->bufs);
    free(c->threads);
    pthread_mutex_destroy(&c->mu);
    pthread_cond_destroy(&c->cv_work);
    pthread_cond_destroy(&c->cv_free);
}

/* Reads exactly need bytes into b. The buffer grows as data actually
 * arrives, so an absurd frame size in a damaged header can't cause an
 * absurd allocation. 0 = ok, 1 = the stream ended early (*got says how
 * far it got), 2 = out of memory. */
static int split_fill(FILE *in, split_buf_t *b, uint64_t need, uint64_t *got_out) {
    size_t got = 0;
    while (got < need) {
        size_t want = (need - got > INGEST_BUF) ? INGEST_BUF : (size_t)(need - got);
        if (b->cap < got + want) {
            size_t ncap = b->cap ? b->cap * 2 : (size_t)INGEST_BUF;
            if (ncap < got + want) ncap = got + want;
            if (ncap > need) ncap = (size_t)need;
            unsigned char *nd = realloc(b->data, ncap);
            if (!nd) { *got_out = got; return 2; }
            b->data = nd;
            b->cap = ncap;
        }
        size_t r = fread(b->data + got, 1, want, in);
        if (r == 0) { *got_out = got; return 1; }
        got += r;
    }
    *got_out = got;
    return 0;
}

typedef struct {
    int alt;                 /* 1 = additive averages, 0/nilsimsa decided by 'nilsimsa' below */
    int nilsimsa;             /* 1 = -n (nilsimsa sort). If neither this nor alt is set, no
                                * sort runs at all: natural filename order, no hashing. */
    int max_sweeps, nthreads, quiet;
    const char *out_path;    /* the one combined sort file; NULL = stdout */
    FILE *dump;               /* -d: already-open handle for additive-average keys
                                * (-a only); opened once by main(), never re-opened.
                                * NULL = don't dump. */
} elem_opts_t;

/* Writes an already-gathered list (mem, memlen) to stdout or to out_path.
 * Only called once every sort has succeeded, so a failed run never
 * touches an existing -f file. Returns 0 or 1. */
static int emit_list(const char *mem, size_t memlen, const char *out_path, size_t lines, int quiet) {
    FILE *out = stdout;
    if (out_path) {
        out = fopen(out_path, "w");
        if (!out) {
            fprintf(stderr, "y4mstore: cannot open %s for writing: %s\n", out_path, strerror(errno));
            return 1;
        }
    }
    int werr = (fwrite(mem, 1, memlen, out) != memlen) || (fflush(out) != 0);
    if (out_path && fclose(out) != 0) werr = 1;
    if (werr) { fprintf(stderr, "y4mstore: write error on the list: %s\n", strerror(errno)); return 1; }
    if (out_path && !quiet) fprintf(stderr, "y4mstore: wrote %zu file path(s) to %s\n", lines, out_path);
    return 0;
}

/* Default (no -n, no -a): no sort runs at all. Exact-duplicate detection
 * is a byproduct of the nilsimsa hashing pass (see nilsimsa_sort_paths
 * above), so it's paired with -n and doesn't run here either -- this
 * just lists the files in natural filename order, with no per-file
 * reads beyond the directory scan itself. Always succeeds. */
static int scan_only(char **paths, size_t count, FILE *out) {
    if (count > 0) {
        qsort(paths, count, sizeof(char *), name_natural_cmp);
        for (size_t i = 0; i < count; i++) fprintf(out, "%s\n", paths[i]);
    }
    return 0;
}

/* Runs the chosen sorter on the y, u and v file lists, in that order, and
 * appends the three sorted lists into one, written to stdout or to
 * op->out_path. The lists are gathered in memory first, so nothing is
 * written (and an existing -f file is left alone) unless all three sorts
 * succeed. sp is only needed, and may be NULL otherwise, for the additive
 * sorter. Returns 0 or 1. */
static int sort_component_lists(char **plist[3], const size_t count[3], const specs_t *sp,
                                const char *specs_dir, const elem_opts_t *op) {
    static const char PLANE[3] = {'y', 'u', 'v'};
    uint64_t psz[3] = {0, 0, 0};
    if (op->alt) plane_sizes(sp, psz);

    char *mem = NULL;
    size_t memlen = 0;
    FILE *ms = open_memstream(&mem, &memlen);
    if (!ms) { fprintf(stderr, "y4mstore: out of memory\n"); return 1; }
    int rc = 1;
    size_t total = 0;
    for (int p = 0; p < 3; p++) {
        if (count[p] == 0) continue;
        char label[2] = {PLANE[p], 0};
        if (!op->quiet) fprintf(stderr, "y4mstore: %s the %c planes (%zu files)\n",
                                op->alt ? "additive-sorting" : op->nilsimsa ? "nilsimsa-sorting" : "listing (no sort)",
                                PLANE[p], count[p]);
        int sr = op->alt
            ? alt_sort(plist[p], count[p], sp, specs_dir, NULL, ms, psz[p], label, op->nthreads, op->quiet, op->dump)
            : op->nilsimsa
                ? nilsimsa_sort_paths(plist[p], count[p], ms, op->max_sweeps, op->nthreads, op->quiet)
                : scan_only(plist[p], count[p], ms);
        if (sr != 0) {
            fprintf(stderr, "y4mstore: sorting the %c planes failed; nothing was written to the list\n", PLANE[p]);
            goto done;
        }
        total += count[p];
    }
    if (fclose(ms) != 0) { ms = NULL; fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
    ms = NULL;

    if (emit_list(mem, memlen, op->out_path, total, op->quiet) != 0) goto done;
    rc = 0;

done:
    if (ms) fclose(ms);
    free(mem);
    return rc;
}

static int elementary_split(const char *src, const char *dest, const elem_opts_t *op) {
    static const char PLANE[3] = {'y', 'u', 'v'};
    const int quiet = op->quiet;
    int rc = 1;
    FILE *in = NULL;
    int pool_started = 0;
    split_ctx_t ctx;
    uint64_t frames = 0;
    char plane_dir[3][PATH_MAX + 64];
    int made_dir[3] = {0, 0, 0}, made_target = 0;
    char named[PATH_MAX];                    /* the subdirectory as given, trailing slashes stripped */
    char target[PATH_MAX + 32];              /* where the planes really go */
    char **plist[3] = {NULL, NULL, NULL};

    if (!strcmp(dest, "-")) {
        fprintf(stderr, "y4mstore: -e needs a subdirectory to write into\n");
        return 1;
    }
    if (strlen(dest) >= PATH_MAX - 64) {
        fprintf(stderr, "y4mstore: the subdirectory path is too long\n");
        return 1;
    }
    off_t fsize = 0;
    if (y4m_open_source(src, "-e", dest, &in, &fsize) != 0) return 1;

    /* 1. the header, before anything is touched */
    specs_t hs;
    uint64_t consumed = 0;
    if (y4m_read_header(in, &hs, quiet, &consumed) != 0) goto done;
    uint64_t psz[3];
    plane_sizes(&hs, psz);
    const uint64_t ysz = psz[0], csz = psz[1];

    /* 2. where the planes go: the named subdirectory, unless it already
     *    holds whole-frame .yuv files -- then <name>-elementary */
    snprintf(named, sizeof(named), "%s", dest);
    for (size_t l = strlen(named); l > 1 && named[l - 1] == '/'; ) named[--l] = 0;
    snprintf(target, sizeof(target), "%s", named);
    {
        struct stat st;
        if (stat(named, &st) == 0 && S_ISDIR(st.st_mode) && has_yuv_files(named)) {
            snprintf(target, sizeof(target), "%s-elementary", named);
            /* shown even with -q: the planes are not where the name given says */
            fprintf(stderr, "y4mstore: %s already holds .yuv frames, so the planes go to %s instead\n", named, target);
        } else if (!quiet) {
            fprintf(stderr, "y4mstore: writing planes to %s\n", target);
        }
    }
    for (int p = 0; p < 3; p++) snprintf(plane_dir[p], sizeof(plane_dir[p]), "%s/%c", target, PLANE[p]);

    /* 3. the y, u and v directories must each be new or empty. A target that
     *    already holds all three planes is the case for sorting, not splitting. */
    if (is_component_dir(target)) {
        fprintf(stderr, "y4mstore: %s already holds y, u and v planes. To sort them, run y4mstore on that directory "
                        "without -e (add -a for the additive sorter); nothing was written.\n", target);
        goto done;
    }
    int to_create[3] = {0, 0, 0};
    for (int p = 0; p < 3; p++) {
        if (dest_new_or_empty(plane_dir[p], &to_create[p],
                              "-e writes only into new or empty plane directories") != 0) goto done;
    }

    /* 4. specs.txt in the set directory (the parent of the subdirectory) */
    char parent[PATH_MAX + 8];
    parent_of(target, parent, sizeof(parent));
    if (specs_settle(parent, &hs, quiet) != 0) goto done;

    /* 5. the directories */
    if (!path_exists(target)) {
        if (mkdir(target, 0777) != 0) { fprintf(stderr, "y4mstore: cannot create %s: %s\n", target, strerror(errno)); goto done; }
        made_target = 1;
    }
    for (int p = 0; p < 3; p++) {
        if (to_create[p]) {
            if (mkdir(plane_dir[p], 0777) != 0 && errno != EEXIST) {
                fprintf(stderr, "y4mstore: cannot create %s: %s\n", plane_dir[p], strerror(errno));
                goto done;
            }
            made_dir[p] = 1;
        }
    }

    /* 6. the frames. This thread reads the stream in order; each whole frame
     *    goes to a worker (-t of them) that writes its y, u and v files. */
    const int nt = op->nthreads < 1 ? 1 : op->nthreads;
    size_t nbufs = (size_t)nt + 2;                                   /* frames in flight */
    {
        uint64_t fit = (256ULL << 20) / (hs.frame_bytes ? hs.frame_bytes : 1);   /* ~256 MiB of buffers at most */
        if (fit < 2) fit = 2;
        if (nbufs > fit) nbufs = (size_t)fit;
    }
    if (split_pool_start(&ctx, nt, nbufs, plane_dir, PLANE, psz) != 0) {
        fprintf(stderr, "y4mstore: cannot start the worker threads\n");
        if (ctx.queue || ctx.threads) split_pool_finish(&ctx);
        goto done;
    }
    pool_started = 1;

    char line[INGEST_LINE_MAX];
    counter_start(g_progress_count && !quiet && fsize == 0);   /* a file source has the progress bar */
    if (g_ctr.on) counter_draw("splitting", 0, 0);              /* visible before the first frame arrives */
    int problem = 0;            /* 0 none, 1 damaged frame header, 2 bad marker, 3 truncated frame, 4 out of memory */
    int problem_fr = 0;
    char problem_shown[24] = "";
    uint64_t problem_got = 0;
    for (;;) {
        int fr = ingest_read_line(in, line, sizeof(line));
        if (fr == 1) break;   /* clean end of stream, on a frame boundary */
        if (fr != 0) { problem = 1; problem_fr = fr; break; }
        if (strncmp(line, "FRAME", 5) != 0 || (line[5] != 0 && line[5] != ' ')) {
            problem = 2;
            safe_text(line, problem_shown, sizeof(problem_shown));
            break;
        }
        consumed += (uint64_t)strlen(line) + 1;

        pthread_mutex_lock(&ctx.mu);                      /* a free buffer, or stop if a worker failed */
        while (ctx.nfree == 0 && !ctx.failed) pthread_cond_wait(&ctx.cv_free, &ctx.mu);
        if (ctx.failed) { pthread_mutex_unlock(&ctx.mu); break; }
        split_buf_t *fb = ctx.freelist[--ctx.nfree];
        pthread_mutex_unlock(&ctx.mu);

        int rr = split_fill(in, fb, hs.frame_bytes, &problem_got);
        if (rr != 0) {                                    /* never hand over a partial frame */
            pthread_mutex_lock(&ctx.mu);
            ctx.freelist[ctx.nfree++] = fb;
            pthread_mutex_unlock(&ctx.mu);
            problem = (rr == 1) ? 3 : 4;
            break;
        }
        pthread_mutex_lock(&ctx.mu);
        ctx.queue[(ctx.qhead + ctx.qcount) % ctx.qcap].buf = fb;
        ctx.queue[(ctx.qhead + ctx.qcount) % ctx.qcap].index = frames + 1;
        ctx.qcount++;
        pthread_cond_signal(&ctx.cv_work);
        pthread_mutex_unlock(&ctx.mu);
        consumed += hs.frame_bytes;
        frames++;
        if (fsize > 0 && !quiet) progress_bar("splitting", (size_t)consumed, (size_t)fsize);
        counter_update("splitting", frames, consumed);
    }
    counter_finish("splitting", frames, consumed);   /* before any message below */
    split_pool_finish(&ctx);   /* every frame handed over is written before we go on */
    pool_started = 0;
    if (ctx.failed) {
        fprintf(stderr, "y4mstore: %s\n", ctx.errmsg);
        goto done;
    }
    if (problem == 1) {
        fprintf(stderr, "y4mstore: stream damaged after %llu frame%s: %s\n", (unsigned long long)frames,
                frames == 1 ? "" : "s", problem_fr == -2 ? "frame header too long" : "input ended inside a frame header");
        goto done;
    }
    if (problem == 2) {
        fprintf(stderr, "y4mstore: stream damaged after %llu frame%s: expected FRAME, found '%s'\n",
                (unsigned long long)frames, frames == 1 ? "" : "s", problem_shown);
        goto done;
    }
    if (problem == 3) {
        char pc = 'y';
        uint64_t in_plane = problem_got;
        if (problem_got >= ysz + csz) { pc = 'v'; in_plane = problem_got - ysz - csz; }
        else if (problem_got >= ysz) { pc = 'u'; in_plane = problem_got - ysz; }
        fprintf(stderr, "y4mstore: stream ended inside frame %llu (%c plane, %llu of %llu bytes); "
                        "%llu complete frame%s kept in %s\n", (unsigned long long)(frames + 1), pc,
                (unsigned long long)in_plane, (unsigned long long)(pc == 'y' ? ysz : csz),
                (unsigned long long)frames, frames == 1 ? "" : "s", target);
        goto done;
    }
    if (problem == 4) { fprintf(stderr, "y4mstore: out of memory reading frame %llu\n", (unsigned long long)(frames + 1)); goto done; }
    if (fsize > 0 && !quiet) progress_bar("splitting", (size_t)fsize, (size_t)fsize);
    if (frames == 0) {
        fprintf(stderr, "y4mstore: the stream has a header but no frames\n");
        goto done;
    }
    if (!quiet) {
        fprintf(stderr, "y4mstore: split %llu frame%s into %s/{y,u,v} (planes of %llu, %llu, %llu bytes)\n",
                (unsigned long long)frames, frames == 1 ? "" : "s", target,
                (unsigned long long)ysz, (unsigned long long)csz, (unsigned long long)csz);
    }

    /* 7. the sort: one plane directory at a time, appended into one list.
     *    The lists are built in frame-number order (not directory order),
     *    so the run is deterministic. */
    size_t counts[3] = {(size_t)frames, (size_t)frames, (size_t)frames};
    for (int p = 0; p < 3; p++) {
        plist[p] = calloc((size_t)frames, sizeof(char *));
        if (!plist[p]) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
        for (uint64_t i = 0; i < frames; i++) {
            char nm[PATH_MAX + 128];
            snprintf(nm, sizeof(nm), "%.*s/yuv-%07llu.%c", PATH_MAX + 63, plane_dir[p], (unsigned long long)(i + 1), PLANE[p]);
            plist[p][i] = strdup(nm);
            if (!plist[p][i]) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
        }
    }
    if (sort_component_lists(plist, counts, &hs, target, op) != 0) {
        fprintf(stderr, "y4mstore: the planes are written in %s but no list was produced\n", target);
        goto done;
    }
    rc = 0;

done:
    counter_break();
    /* a directory this run created but never filled is not left behind */
    if (rc != 0 && frames == 0) {
        for (int p = 2; p >= 0; p--) if (made_dir[p]) rmdir(plane_dir[p]);
        if (made_target) rmdir(target);
    }
    for (int p = 0; p < 3; p++) {
        if (plist[p]) {
            for (uint64_t i = 0; i < frames; i++) free(plist[p][i]);
            free(plist[p]);
        }
    }
    (void)pool_started;
    if (in && in != stdin) fclose(in);
    return rc;
}

/* ---------------------------------------------------------------- */
/* Sorting an existing 3-component directory (no -e)                 */
/*                                                                   */
/* A directory that already holds y/, u/ and v/ planes (what -e      */
/* writes) is recognised on its own: the selected sorter runs on each */
/* component directory and the three lists are appended into one, y   */
/* then u then v -- the same sort phase -e ends with, without the     */
/* split. The specs come from the parent (set) directory, like a      */
/* clip's.                                                            */
/* ---------------------------------------------------------------- */

static int component_sort(const char *dir, const elem_opts_t *op) {
    static const char PLANE[3] = {'y', 'u', 'v'};
    static const char *const EXT[3] = {".y", ".u", ".v"};
    int rc = 1;
    char **plist[3] = {NULL, NULL, NULL};
    size_t count[3] = {0, 0, 0};

    if (!op->quiet) {
        fprintf(stderr, "y4mstore: %s -- component storage (y, u, v planes); sort: %s\n", dir,
                op->alt ? "additive averages (-a)" : op->nilsimsa ? "nilsimsa (-n)" : "none -- natural order list being generated");
    }
    /* the specs are found the way a clip's are (its parent's specs.txt). The
     * additive sorter needs them for the plane sizes; nilsimsa does not, so
     * it only reads one that is already there. */
    char sdir[PATH_MAX + 8];
    specs_locate(dir, sdir, sizeof(sdir), op->quiet);
    char spath[PATH_MAX + 32];
    snprintf(spath, sizeof(spath), "%s/%s", sdir, SPECS_FILENAME);
    specs_t sp;
    memset(&sp, 0, sizeof(sp));
    int have_specs = 0;
    if (path_exists(spath) || op->alt) {
        if (specs_for_dir(sdir, &sp, op->quiet) != 0) return 1;
        have_specs = 1;
    }

    for (int p = 0; p < 3; p++) {
        char sub[PATH_MAX + 8];
        snprintf(sub, sizeof(sub), "%.*s/%c", PATH_MAX - 4, dir, PLANE[p]);
        plist[p] = list_ext_files(sub, EXT[p], &count[p]);
        if (!plist[p] || count[p] == 0) {
            fprintf(stderr, "y4mstore: no %s files found in %s\n", EXT[p], sub);
            goto done;
        }
    }
    if ((count[0] != count[1] || count[1] != count[2]) && !op->quiet) {
        fprintf(stderr, "y4mstore: warning: the component directories hold different numbers of files "
                        "(y=%zu, u=%zu, v=%zu); each is sorted as it is\n", count[0], count[1], count[2]);
    }
    rc = sort_component_lists(plist, count, have_specs ? &sp : NULL, sdir, op);

done:
    for (int p = 0; p < 3; p++) {
        if (plist[p]) {
            for (size_t i = 0; i < count[p]; i++) free(plist[p][i]);
            free(plist[p]);
        }
    }
    return rc;
}

/* ---------------------------------------------------------------- */
/* Several directories, one sort file (e.g. multiple episodes)       */
/*                                                                   */
/*   y4mstore [-a] DIR DIR ...                                        */
/*                                                                   */
/* For sorting only. Every directory must be the same kind: all      */
/* whole-frame clips, or all y/u/v component directories. Frame case: */
/* all the frame files are pooled and sorted once into one list.      */
/* Component case: all the y planes are pooled and sorted, then all   */
/* the u planes, then all the v planes, and the three lists appended. */
/* Directories are taken in the order given, so a glob (ep*) gives    */
/* episode order.                                                     */
/* ---------------------------------------------------------------- */

typedef enum { DK_OTHER = 0, DK_FRAMES, DK_COMPONENT } dirkind_t;

static dirkind_t classify_dir(const char *dir) {
    if (is_component_dir(dir)) return DK_COMPONENT;
    if (count_stray_top_level_files(dir) > 0) return DK_FRAMES;
    return DK_OTHER;
}

static const char *dirkind_name(dirkind_t k) {
    switch (k) {
    case DK_FRAMES: return "whole-frame clip";
    case DK_COMPONENT: return "y/u/v component directory";
    default: return "neither (no frame files and no y/u/v planes)";
    }
}

typedef struct {
    const char *dir;
    dirkind_t kind;
    dev_t dev;
    ino_t ino;
    char sdir[PATH_MAX + 8];   /* the directory whose specs.txt applies */
    specs_t specs;
    int have_specs;
    int specs_done;
    size_t nfiles;             /* frame files (frame clips) */
    size_t nplane[3];          /* y, u, v files (component directories) */
} multi_dir_t;

/* Only the fields that change what a sample means: additive keys from an
 * 8-bit and a 10-bit stream, or from different geometries, aren't comparable. */
static int geom_equal(const specs_t *a, const specs_t *b) {
    return a->width == b->width && a->height == b->height && a->sampling == b->sampling && a->bitdepth == b->bitdepth;
}

static void geom_text(const specs_t *sp, char *buf, size_t n) {
    snprintf(buf, n, "%dx%d %s %d-bit", sp->width, sp->height, SAMPLINGS[sp->sampling].name, sp->bitdepth);
}

/* Appends the k strings of src to *dst (taking ownership). 0 or -1. */
static int append_paths(char ***dst, size_t *n, size_t *cap, char **src, size_t k) {
    if (*n + k > *cap) {
        size_t nc = *cap ? *cap : 1024;
        while (nc < *n + k) nc *= 2;
        char **t = realloc(*dst, nc * sizeof(char *));
        if (!t) return -1;
        *dst = t;
        *cap = nc;
    }
    memcpy(*dst + *n, src, k * sizeof(char *));
    *n += k;
    return 0;
}

static int multi_sort(int ndirs, char **dirs, const elem_opts_t *op) {
    static const char PLANE[3] = {'y', 'u', 'v'};
    static const char *const EXT[3] = {".y", ".u", ".v"};
    int rc = 1;
    const int quiet = op->quiet;
    char **all = NULL;                          /* frame case: every frame file */
    size_t nall = 0, capall = 0;
    char **pl[3] = {NULL, NULL, NULL};          /* component case: every y / u / v file */
    size_t np[3] = {0, 0, 0}, cp[3] = {0, 0, 0};
    multi_dir_t *md = calloc((size_t)ndirs, sizeof(*md));
    if (!md) { fprintf(stderr, "y4mstore: out of memory\n"); return 1; }

    /* 1. each must be an existing directory, and none may be given twice */
    for (int i = 0; i < ndirs; i++) {
        struct stat st;
        md[i].dir = dirs[i];
        if (stat(dirs[i], &st) != 0) { fprintf(stderr, "y4mstore: cannot stat %s: %s\n", dirs[i], strerror(errno)); goto done; }
        if (!S_ISDIR(st.st_mode)) { fprintf(stderr, "y4mstore: %s is not a directory\n", dirs[i]); goto done; }
        md[i].dev = st.st_dev;
        md[i].ino = st.st_ino;
        for (int j = 0; j < i; j++) {
            if (md[j].dev == md[i].dev && md[j].ino == md[i].ino) {
                fprintf(stderr, "y4mstore: %s and %s are the same directory; give each directory once "
                                "(its files would be listed twice)\n", dirs[j], dirs[i]);
                goto done;
            }
        }
    }

    /* 2. the sanity check: they must all be the same kind of directory */
    int bad = 0, mixed = 0;
    for (int i = 0; i < ndirs; i++) {
        md[i].kind = classify_dir(dirs[i]);
        if (md[i].kind == DK_OTHER) bad = 1;
        if (md[i].kind != md[0].kind) mixed = 1;
    }
    if (bad || mixed) {
        fprintf(stderr, "y4mstore: %s. Several directories can be sorted together only if every one is a "
                        "whole-frame clip, or every one is a y/u/v component directory:\n",
                bad ? "some directories are neither" : "the directories are not all the same type");
        for (int i = 0; i < ndirs; i++) fprintf(stderr, "  %-40s %s\n", dirs[i], dirkind_name(md[i].kind));
        fprintf(stderr, "         Nothing was written.\n");
        goto done;
    }
    const dirkind_t kind = md[0].kind;

    /* 3. specs, found the way a single directory's are (the parent's specs.txt).
     *    Directories in one set share theirs, so each file is read once. The
     *    additive sorter needs them for every directory; nilsimsa reads them for
     *    frame clips (as it does for one directory) and otherwise only if present. */
    for (int i = 0; i < ndirs; i++) {
        specs_locate(md[i].dir, md[i].sdir, sizeof(md[i].sdir), quiet);
        for (int j = 0; j < i; j++) {
            if (md[j].specs_done && !strcmp(md[j].sdir, md[i].sdir)) {
                md[i].specs = md[j].specs;
                md[i].have_specs = md[j].have_specs;
                md[i].specs_done = 1;
                break;
            }
        }
        if (md[i].specs_done) continue;
        char spath[PATH_MAX + 32];
        snprintf(spath, sizeof(spath), "%s/%s", md[i].sdir, SPECS_FILENAME);
        if (kind == DK_FRAMES || op->alt || path_exists(spath)) {
            if (specs_for_dir(md[i].sdir, &md[i].specs, quiet) != 0) goto done;
            md[i].have_specs = 1;
        }
        md[i].specs_done = 1;
    }
    int first = -1;
    for (int i = 0; i < ndirs; i++) {
        if (!md[i].have_specs) continue;
        if (first < 0) { first = i; continue; }
        if (!geom_equal(&md[first].specs, &md[i].specs)) {
            char a[64], b[64];
            geom_text(&md[first].specs, a, sizeof(a));
            geom_text(&md[i].specs, b, sizeof(b));
            if (op->alt) {
                fprintf(stderr, "y4mstore: the directories do not share one format, so their additive keys are not "
                                "comparable:\n  %-40s %s\n  %-40s %s\n         Nothing was written.\n",
                        md[first].dir, a, md[i].dir, b);
                goto done;
            }
            if (!quiet) {
                fprintf(stderr, "y4mstore: warning: %s is %s but %s is %s; they are sorted together anyway\n",
                        md[first].dir, a, md[i].dir, b);
            }
        }
    }

    /* 4. walk each directory, in the order given, pooling its files */
    for (int i = 0; i < ndirs; i++) {
        if (kind == DK_FRAMES) {
            size_t c = 0;
            char **l = list_ext_files(md[i].dir, NULL, &c);
            if (!l || c == 0) { fprintf(stderr, "y4mstore: no frame files found in %s\n", md[i].dir); free(l); goto done; }
            if (append_paths(&all, &nall, &capall, l, c) != 0) {
                for (size_t k = 0; k < c; k++) free(l[k]);
                free(l);
                fprintf(stderr, "y4mstore: out of memory\n");
                goto done;
            }
            free(l);
            md[i].nfiles = c;
        } else {
            for (int p = 0; p < 3; p++) {
                char sub[PATH_MAX + 8];
                snprintf(sub, sizeof(sub), "%.*s/%c", PATH_MAX - 4, md[i].dir, PLANE[p]);
                size_t c = 0;
                char **l = list_ext_files(sub, EXT[p], &c);
                if (!l || c == 0) { fprintf(stderr, "y4mstore: no %s files found in %s\n", EXT[p], sub); free(l); goto done; }
                if (append_paths(&pl[p], &np[p], &cp[p], l, c) != 0) {
                    for (size_t k = 0; k < c; k++) free(l[k]);
                    free(l);
                    fprintf(stderr, "y4mstore: out of memory\n");
                    goto done;
                }
                free(l);
                md[i].nplane[p] = c;
            }
            if (!quiet && (md[i].nplane[0] != md[i].nplane[1] || md[i].nplane[1] != md[i].nplane[2])) {
                fprintf(stderr, "y4mstore: warning: %s holds different numbers of plane files (y=%zu, u=%zu, v=%zu)\n",
                        md[i].dir, md[i].nplane[0], md[i].nplane[1], md[i].nplane[2]);
            }
        }
    }

    const char *sortname = op->alt ? "additive averages (-a)" : op->nilsimsa ? "nilsimsa (-n)" : "none -- natural order list being generated";
    if (!quiet) {
        if (kind == DK_FRAMES) {
            fprintf(stderr, "y4mstore: %d directories -- frame storage; sort: %s; %zu files pooled into one list\n",
                    ndirs, sortname, nall);
            for (int i = 0; i < ndirs; i++) fprintf(stderr, "  %-40s %zu files\n", md[i].dir, md[i].nfiles);
        } else {
            fprintf(stderr, "y4mstore: %d directories -- component storage (y, u, v planes); sort: %s; "
                            "pooling the y planes, then the u, then the v (%zu, %zu, %zu files)\n",
                    ndirs, sortname, np[0], np[1], np[2]);
            for (int i = 0; i < ndirs; i++)
                fprintf(stderr, "  %-40s y=%zu u=%zu v=%zu\n", md[i].dir, md[i].nplane[0], md[i].nplane[1], md[i].nplane[2]);
        }
    }

    /* 5. the sort: one pooled sort per list, appended into one sort file */
    const specs_t *sp = first >= 0 ? &md[first].specs : NULL;
    const char *specs_dir = first >= 0 ? md[first].sdir : md[0].sdir;
    if (kind == DK_COMPONENT) {
        rc = sort_component_lists(pl, np, sp, specs_dir, op);
    } else if (op->alt) {
        rc = alt_sort(all, nall, sp, specs_dir, op->out_path, NULL, 0, NULL, op->nthreads, quiet, op->dump);
    } else {
        char *mem = NULL;
        size_t memlen = 0;
        FILE *ms = open_memstream(&mem, &memlen);
        if (!ms) { fprintf(stderr, "y4mstore: out of memory\n"); goto done; }
        int sr = op->nilsimsa
            ? nilsimsa_sort_paths(all, nall, ms, op->max_sweeps, op->nthreads, quiet)
            : scan_only(all, nall, ms);
        if (fclose(ms) != 0) sr = 1;
        rc = sr == 0 ? emit_list(mem, memlen, op->out_path, nall, quiet) : 1;
        free(mem);
    }

done:
    for (size_t i = 0; i < nall; i++) free(all[i]);
    free(all);
    for (int p = 0; p < 3; p++) {
        for (size_t i = 0; i < np[p]; i++) free(pl[p][i]);
        free(pl[p]);
    }
    free(md);
    return rc;
}

/* ---------------------------------------------------------------- */
/* main                                                              */
/* ---------------------------------------------------------------- */

static void usage(FILE *fp, const char *prog) {
    fprintf(fp,
        "Usage: %s [options] <dir> [<dir> ...]\n"
        "       find /some/dir -type f | %s [options] -\n"
        "       %s -s [-f FILE] <clipdir>      (stream a clip as yuv4mpeg)\n"
        "       %s -c <setdir>                 (check/create specs.txt)\n"
        "       %s -n [-t N] [-f FILE] <dir>   (nilsimsa similarity sort)\n"
        "       %s -a [-t N] [-f FILE] [-d FILE] <dir>   (sort by Y/U/V sample averages)\n"
        "       %s -i FILE|- <clipdir>         (read a y4m stream into a clip)\n"
        "       %s -e FILE|- [-n|-a] [-f FILE] <subdir>  (split a y4m stream into y/u/v planes, sort each)\n\n"
        "<dir> holds specs.txt (resolution, sampling, framerate, bitdepth, aspect) and\n"
        "one subdirectory per clip; if specs.txt is missing you are asked for it.\n"
        "A <dir> that already holds y/, u/ and v/ planes (made by -e) is sorted per\n"
        "component, without -e: each plane directory is sorted and the lists appended.\n"
        "Several directories (say, several episodes) sort together into one list, but only\n"
        "if they are all whole-frame clips or all component directories; component\n"
        "directories pool all the y planes, then the u, then the v.\n\n"
        "A <dir> that holds only subdirectories (a set with more than one clip in it) is\n"
        "refused rather than blindly recursed into every clip as one undifferentiated\n"
        "pile: name or glob the clips explicitly (see above) instead.\n\n"
        "Without -n or -a, no sort runs: files are listed in natural filename order,\n"
        "with no hashing and no duplicate detection (that's paired with -n, see below).\n"
        "This is the default so that adding clips over time, or checking a growing\n"
        "pooled set, doesn't force a full sort every run -- pass -n when you're ready\n"
        "to actually sort and archive. Before any of that starts, a line is printed\n"
        "naming the storage layout (frame or component) and the chosen sort (or\n"
        "\"none\"), so you can Ctrl+C if it's not what you meant to run.\n\n"
        "-f and -d never append: if the file already exists, you're asked whether to\n"
        "overwrite it (refused outright if stdin isn't a terminal, rather than guessing).\n\n"
        "Options:\n"
        "  -n     nilsimsa similarity sort (this used to be the unconditional\n"
        "         default). Exact-duplicate detection (a byte-identical-file\n"
        "         check) is a byproduct of this hashing pass, so it only runs\n"
        "         with -n. Mutually exclusive with -a.\n"
        "  -o N   max nilsimsa passes (default 50); -n only\n"
        "  -t N   worker threads, per frame: hashing (-n), averaging\n"
        "         (-a), and writing the y/u/v plane files (-e). Default = number\n"
        "         of CPUs\n"
        "  -i SRC ingest: read a yuv4mpeg (y4m) stream from SRC (a file, or - for\n"
        "         stdin) and write each frame as a raw file (yiv-000001.yuv, ...)\n"
        "         into the new or empty clip directory <clipdir>. specs.txt in its\n"
        "         parent is created from the stream header; if one exists it must\n"
        "         match, otherwise nothing is written.\n"
        "  -e SRC elementary: read a y4m stream from SRC (a file, or - for stdin;\n"
        "         not a directory -- pipe another instance in: y4mstore -s DIR |\n"
        "         y4mstore -e - NAME). Each frame is cut into planes written as\n"
        "         <subdir>/y/yuv-0000001.y, u/...u, v/...v (or <subdir>-elementary\n"
        "         if <subdir> already holds .yuv files). Then the chosen sorter\n"
        "         (nilsimsa with -n, additive averages with -a, or, by default,\n"
        "         no sort at all) runs on the y, u and v planes in that order,\n"
        "         and the three lists are appended into one (-f FILE\n"
        "         or stdout).\n"
        "  -a     alternate sort, no similarity hash: each frame's Y, U and V\n"
        "         sample averages (scaled x100 for precision) become a key,\n"
        "         Y first, most significant, digit width per bit depth\n"
        "         (5 digits/plane at 8-bit, 6 at 10-bit). Needs specs.txt and\n"
        "         exactly one frame per file. -t applies; -o does not.\n"
        "         Equal keys are ordered by filename.\n"
        "  -d PATH with -a: outputs averaged per frame/component to PATH\n"
        "         (sorted, one per line, as `KEY  path`). Requires -a.\n"
        "  -c     check: look for specs.txt in <setdir>. If it is missing, ask\n"
        "         for the values and create it; if present, list the specs.\n"
        "         Nothing is scanned, sorted or served.\n"
        "  -s     serve: write <clipdir> as a yuv4mpeg (y4m) stream on stdout,\n"
        "         (a directory of y/, u/ and v/ planes from -e is served back the\n"
        "         same way, each frame's planes paired by filename),\n"
        "         using specs.txt from its parent directory. Pipe it to ffplay,\n"
        "         mpv or ffmpeg, or redirect to a .y4m file. Frames go out in\n"
        "         filename order (numbers compare as numbers). -o/-t\n"
        "         do not apply.\n"
        "  -f PATH write the file list (or, with -s, the y4m stream)\n"
        "         to PATH instead of stdout\n"
        "  -p     with -i or -e reading stdin (-): show a running count of frames\n"
        "         read, so a slow upstream doesn't look like a hang (a file source\n"
        "         already has a progress bar). One redrawn line on a terminal, a\n"
        "         plain line every few seconds if stderr is redirected. -q wins.\n"
        "  -q     quiet (suppress progress messages on stderr)\n"
        "  -h     show this help\n",
        prog, prog, prog, prog, prog, prog, prog, prog);
}

int main(int argc, char **argv) {
    int max_sweeps = 50;
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1) nthreads = 1;
    int quiet = 0;
    int serve = 0;
    int check = 0;
    int alt = 0;
    int do_nilsimsa = 0;
    const char *ingest_src = NULL;
    const char *elem_src = NULL;
    const char *out_path = NULL;
    const char *dump_path = NULL;
    int opt;

    while ((opt = getopt(argc, argv, "o:t:f:i:e:d:acsn" "pqh")) != -1) {
        switch (opt) {
            case 'o': max_sweeps = atoi(optarg); break;
            case 't': nthreads = atoi(optarg); break;
            case 'f': out_path = optarg; break;
            case 'd': dump_path = optarg; break;
            case 's': serve = 1; break;
            case 'c': check = 1; break;
            case 'a': alt = 1; break;
            case 'n': do_nilsimsa = 1; break;
            case 'i': ingest_src = optarg; break;
            case 'e': elem_src = optarg; break;
            case 'p': g_progress_count = 1; break;
            case 'q': quiet = 1; break;
            case 'h': usage(stdout, argv[0]); return 0;   /* asked for: stdout, so it can be piped */
            default: usage(stderr, argv[0]); return 1;
        }
    }
    if (optind >= argc) {
        usage(stderr, argv[0]);
        return 1;
    }
    if (alt && do_nilsimsa) {
        fprintf(stderr, "y4mstore: -a and -n are mutually exclusive (pick one sorter)\n");
        return 1;
    }
    if (dump_path && !alt) {
        fprintf(stderr, "y4mstore: -d requires -a (it dumps additive-average keys)\n");
        return 1;
    }

    /* -f and -d are never appended to: if either already exists, ask
     * before truncating it (or refuse outright if not on a terminal).
     * -d is opened exactly once, here, and the same handle is reused for
     * every plane sorted in this run (y, u, v) rather than reopened. */
    if (out_path && !confirm_overwrite(out_path, quiet)) return 1;
    FILE *dump_fp = NULL;
    if (dump_path) {
        if (!confirm_overwrite(dump_path, quiet)) return 1;
        dump_fp = fopen(dump_path, "w");
        if (!dump_fp) {
            fprintf(stderr, "y4mstore: cannot open %s for writing: %s\n", dump_path, strerror(errno));
            return 1;
        }
    }

    /* Several directories are for sorting only. The other modes take exactly
     * one, and must say so rather than quietly ignore the rest. */
    const int ndirs = argc - optind;
    if (ndirs > 1) {
        if (ingest_src || elem_src || serve || check) {
            fprintf(stderr, "y4mstore: %s takes exactly one directory; several directories are for sorting only "
                            "(y4mstore [-a] DIR DIR ...)\n", ingest_src ? "-i" : elem_src ? "-e" : serve ? "-s" : "-c");
            return 1;
        }
        for (int i = optind; i < argc; i++) {
            if (!strcmp(argv[i], "-")) {
                fprintf(stderr, "y4mstore: a file list on stdin (-) can't be combined with directories\n");
                return 1;
            }
        }
    }

    /* Printed unconditionally (even with -q) and before any real work
     * starts, so a forgotten -f is visible right away -- Ctrl+C here
     * costs nothing, discovering it after a long run costs everything. */
    if (ingest_src && (alt || do_nilsimsa || serve || check || elem_src)) {
        fprintf(stderr, "y4mstore: -i cannot be combined with -a, -n, -s, -c or -e\n");
        return 1;
    }
    if (elem_src && (serve || check)) {
        fprintf(stderr, "y4mstore: -e cannot be combined with -s or -c (it can be combined with -a or -n: additive or nilsimsa sorter)\n");
        return 1;
    }
    if (elem_src) {
        fprintf(stderr, "y4mstore: elementary: y4m from %s into %s -- component storage (y, u, v planes); sort: %s\n",
                !strcmp(elem_src, "-") ? "stdin" : elem_src, argv[optind],
                alt ? "additive averages (-a)" : do_nilsimsa ? "nilsimsa (-n)" : "none -- natural order list being generated");
        elem_opts_t eo = { alt, do_nilsimsa, max_sweeps, nthreads, quiet, out_path, dump_fp };
        int erc = elementary_split(elem_src, argv[optind], &eo);
        if (dump_fp) fclose(dump_fp);
        return erc;
    }
    if (alt && (serve || check)) {
        fprintf(stderr, "y4mstore: -a cannot be combined with -s or -c\n");
        return 1;
    }
    if (do_nilsimsa && (serve || check)) {
        fprintf(stderr, "y4mstore: -n cannot be combined with -s or -c\n");
        return 1;
    }
    if (ingest_src) {
        fprintf(stderr, "y4mstore: ingest: y4m from %s into %s\n",
                !strcmp(ingest_src, "-") ? "stdin" : ingest_src, argv[optind]);
        return ingest_y4m(ingest_src, argv[optind], quiet);
    }
    if (check) {
        if (serve) {
            fprintf(stderr, "y4mstore: -c and -s cannot be combined\n");
            return 1;
        }
        if (!strcmp(argv[optind], "-")) {
            fprintf(stderr, "y4mstore: -c needs a set directory, not a file list on stdin\n");
            return 1;
        }
        struct stat cst;
        if (stat(argv[optind], &cst) != 0) {
            fprintf(stderr, "y4mstore: cannot stat %s: %s\n", argv[optind], strerror(errno));
            return 1;
        }
        if (!S_ISDIR(cst.st_mode)) {
            fprintf(stderr, "y4mstore: -c needs a set directory; %s is not a directory\n", argv[optind]);
            return 1;
        }
        return check_specs(argv[optind], quiet);
    }

    if (serve) {
        fprintf(stderr, "y4mstore: serve: y4m stream from %s | outputting to %s\n",
                argv[optind], out_path ? out_path : "stdout");
        if (!strcmp(argv[optind], "-")) {
            fprintf(stderr, "y4mstore: -s needs a clip directory, not a file list on stdin\n");
            return 1;
        }
        struct stat sst;
        if (stat(argv[optind], &sst) != 0) {
            fprintf(stderr, "y4mstore: cannot stat %s: %s\n", argv[optind], strerror(errno));
            return 1;
        }
        if (!S_ISDIR(sst.st_mode)) {
            fprintf(stderr, "y4mstore: -s needs a clip directory; %s is not a directory\n", argv[optind]);
            return 1;
        }
        return serve_clip(argv[optind], out_path, quiet);
    }

    /* Several directories (say, several episodes): pooled into one sort file */
    if (ndirs > 1) {
        elem_opts_t mo = { alt, do_nilsimsa, max_sweeps, nthreads, quiet, out_path, dump_fp };
        int mrc = multi_sort(ndirs, &argv[optind], &mo);
        if (dump_fp) fclose(dump_fp);
        return mrc;
    }

    /* No -e given, but the directory already holds y/u/v planes: sort each
     * component directory and append the lists. (Decided before the output
     * file is opened, so nothing is truncated for nothing.) */
    if (strcmp(argv[optind], "-") != 0 && is_component_dir(argv[optind])) {
        elem_opts_t co = { alt, do_nilsimsa, max_sweeps, nthreads, quiet, out_path, dump_fp };
        int crc = component_sort(argv[optind], &co);
        if (dump_fp) fclose(dump_fp);
        return crc;
    }

    /* Not a clip itself, and not a component directory: if it holds only
     * subdirectories (no loose frame files directly inside), it's a set
     * directory with more than one clip in it. Pointing a sort at it
     * directly would otherwise recurse blindly through every clip's files
     * as one undifferentiated pile -- mixing separate clips together, and
     * for component clips, mixing y/u/v planes together too. Refuse
     * rather than silently doing that. */
    if (strcmp(argv[optind], "-") != 0 && dir_holds_only_subdirs(argv[optind])) {
        fprintf(stderr, "y4mstore: multiple clips detected. Please list the clips to be ingested by name or glob\n");
        if (dump_fp) fclose(dump_fp);
        return 1;
    }

    fprintf(stderr, "y4mstore: %s -- frame storage; sort: %s | outputting to %s\n",
            argv[optind],
            alt ? "additive averages (-a)"
                : do_nilsimsa ? "nilsimsa (-n)"
                              : "none -- natural order list being generated (pass -n for nilsimsa or -a for additive)",
            out_path ? out_path : "stdout");

    FILE *out = stdout;
    if (out_path && !alt) {   /* -a opens its output only once every frame is keyed */
        out = fopen(out_path, "w");
        if (!out) {
            fprintf(stderr, "y4mstore: cannot open %s for writing: %s\n",
                    out_path, strerror(errno));
            return 1;
        }
    }

    const char *src = argv[optind];
    pathlist_t pl;
    pathlist_init(&pl);
    char alt_sdir[PATH_MAX + 8] = "";   /* where specs.txt was found (for -a) */
    specs_t alt_specs;

    if (alt && !strcmp(src, "-")) {
        fprintf(stderr, "y4mstore: -a needs a directory (its frame size comes from specs.txt), not a file list on stdin\n");
        return 1;
    }
    if (!strcmp(src, "-")) {
        read_stdin_list(&pl);
    } else {
        struct stat st;
        if (stat(src, &st) != 0) {
            fprintf(stderr, "y4mstore: cannot stat %s: %s\n", src, strerror(errno));
            return 1;
        }
        if (S_ISDIR(st.st_mode)) {
            char sdir[PATH_MAX + 8];   /* the directory whose specs.txt applies */
            specs_locate(src, sdir, sizeof(sdir), quiet);
            specs_t specs;
            if (specs_for_dir(sdir, &specs, quiet) != 0) return 1;
            alt_specs = specs;
            snprintf(alt_sdir, sizeof(alt_sdir), "%s", sdir);
            scan_dir(src, &pl, quiet, 1);
        } else if (alt) {
            fprintf(stderr, "y4mstore: -a needs a directory; %s is not a directory\n", src);
            return 1;
        } else {
            pathlist_push(&pl, src);
        }
    }

    if (pl.count == 0) {
        if (!quiet) fprintf(stderr, "y4mstore: no files found\n");
        return 0;
    }

    if (alt) {
        int arc = alt_sort(pl.paths, pl.count, &alt_specs, alt_sdir, out_path, NULL, 0, NULL, nthreads, quiet, dump_fp);
        for (size_t i = 0; i < pl.count; i++) free(pl.paths[i]);
        free(pl.paths);
        if (dump_fp) fclose(dump_fp);
        return arc;
    }

    int nrc = do_nilsimsa
        ? nilsimsa_sort_paths(pl.paths, pl.count, out, max_sweeps, nthreads, quiet)
        : scan_only(pl.paths, pl.count, out);
    if (nrc == 0 && out_path) {
        fclose(out);
        if (!quiet) fprintf(stderr, "y4mstore: wrote %zu file path(s) to %s\n", pl.count, out_path);
    }
    for (size_t i = 0; i < pl.count; i++) free(pl.paths[i]);
    free(pl.paths);
    return nrc;
}
