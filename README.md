# y4mstore

Stores raw YUV frame data ingested from a y4m stream on disk, in a
set/clip/plane directory layout with per-set `specs.txt`, and serves it
back out as y4m on demand. It can also sort the stored frames by
[nilsimsa](https://en.wikipedia.org/wiki/Nilsimsa_Hash) similarity before
archiving, so that similar frames end up adjacent in the output list,
which can improve compression on redundant data — but storing and
round-tripping y4m is the main job; sorting is one thing you can do with
what's stored.

```
$ ffmpeg -i ep1.mkv -pix_fmt yuv420p -f yuv4mpegpipe - | y4mstore -i - myset/ep1
$ y4mstore -s myset/ep1 | ffplay -
$ y4mstore -n myset/ep1 | tar -T- --no-recursion -cf - | zstd --long -o frames.tar.zst
```

y4mstore is a spinoff of **nilsort**, narrowed to raw (uncompressed) YUV. The MinHash hasher (`-M`) and the EBML/Matroska container handling (`-m`) were removed: nilsort's benchmarks showed plain nilsimsa clearly beating MinHash on raw YUV, and raw frames have no container framing to strip. `-i` turns an incoming y4m stream into a stored clip directory, `-s` serves a stored clip directory back out as `yuv4mpeg` (y4m) so it can be played or converted, and `-e` splits a y4m stream into separate Y, U and V plane files and sorts each set.

## How it works

Sorting is opt-in, via `-n` (nilsimsa) or `-a` (additive averages) — two
different routes to a file order, described separately below. Without
either, no sort runs at all: files are listed in natural filename order,
with no hashing and no duplicate detection — that way, adding clips over
time, or just checking a growing pooled set, doesn't force a full sort on
every run. Pass `-n` or `-a` once, when you're ready to actually sort and
archive.

### The nilsimsa route (`-n`)

1. **Hash.** Every file gets a 256-bit nilsimsa digest (multithreaded — one thread pool hashes all files in parallel, each digest independent of the others). Exact-duplicate detection (see [Duplicate detection](#duplicate-detection)) is a byproduct of this same pass, so it only runs with `-n`.
2. **Order.** Files are greedily chained: start anywhere, repeatedly append whichever unvisited file is most similar to the last one added.
3. **Refine.** [2-opt](https://en.wikipedia.org/wiki/2-opt) local-search passes look for pairs of segments to reverse if doing so increases total adjacent similarity. Runs until it stops finding improvements (or hits a sweep cap).

By default the **entire input is treated as one batch** — every file is weighed against every other for maximum accuracy. Only digests (32 bytes/file) and index arrays are held in memory, so the O(n²) cost is comparison *time*, not memory.

This ordering step is the one meant to be handed to zpaq-order (see [Intended pairing](#intended-pairing-zpaq--zpaqfuse-and-a-caveat-for-hdds) below): `-n` decides the order, zpaq-order compresses the files in exactly that order, and any exact-duplicate work — the byte-level count described in [Duplicate detection](#duplicate-detection), and the real fragment-level dedup that matters for storage — happens separately, downstream, once zpaq-order runs. Nothing about the similarity ordering itself detects or relies on duplicates.

### The additive-average route (`-a`)

A much cheaper alternative: no similarity hash, no O(n²) comparisons, and no per-file duplicate detection.

1. **Split** each frame into its Y, U and V planes, using the resolution, sampling and bit depth from `specs.txt`.
2. **Average** each plane's samples, rounded to the nearest whole number.
3. **Key** each frame as a 12-digit number (4 digits per plane, Y first, most significant). In component mode, each file is only one plane, so only a 4-digit number is assigned per file.
4. **Sort** frames numerically by that key, ascending, with ties broken by filename.

Every frame is handled independently, so this is O(n): a fixed thread pool averages frames in parallel, and the whole run ends in a single numeric sort, with no refinement passes. It's not a similarity measure — two frames with the same three averages can look nothing alike — and it cannot detect exact duplicates, since averaging is lossy by nature. Full details, the caveats, and measured numbers are in [Alternate hash: additive averages](#alternate-hash-additive-averages--a) below.

## Build

```bash
gcc -O2 -Wall -pthread -march=native -o y4mstore y4mstore.c
```

`-march=native` builds for the exact CPU you're compiling on and picks up hardware POPCNT automatically (digest comparison is a Hamming-distance calculation; POPCNT measured ~8x faster than a lookup table in nilsort). Build on each machine you run it on. To build once and run on other x86-64 machines, use `-mpopcnt` instead of `-march=native`.

## Usage

```
y4mstore [options] <dir>
find /some/dir -type f | y4mstore [options] -
y4mstore -s [-f FILE] <clipdir or y/u/v planes dir>
y4mstore -c <setdir>
y4mstore -n [-f FILE] <dir>
y4mstore -a [-f FILE] [-d FILE] <dir>
y4mstore -i FILE|- <clipdir>
y4mstore -e FILE|- [-n|-a] [-f FILE] <subdir>
y4mstore [-n|-a] [-f FILE] <dir with y/u/v planes>     (per-component sort, no -e)
y4mstore [-n|-a] [-f FILE] DIR DIR DIR ...             (several episodes, one sort file)
```

Without `-n` or `-a`, no sort runs: the file list is printed in natural filename order, with no hashing and no duplicate detection (that's paired with `-n`; see [How it works](#how-it-works)). This is the default so that adding clips over time, or checking a growing pooled set, doesn't force a full sort on every run. Before any of that starts, a line names the storage layout (frame or component) and the chosen sort (or "none"), so there's a window to Ctrl+C if it's not what you meant to run.

`-f` and `-d` never append — each run overwrites fully, start to finish. If the target file already exists, you're asked whether to overwrite it; off a terminal (piped, scripted) this is refused outright rather than guessed at.

| Option | Meaning |
|---|---|
| `-n` | **Nilsimsa sort:** the similarity sort described in [How it works](#how-it-works). Exact-duplicate detection is a byproduct of this pass, so it only runs with `-n`. Mutually exclusive with `-a`. |
| `-a` | **Alternate hash:** sort by additive Y/U/V sample averages instead of nilsimsa similarity. See [Alternate hash: additive averages](#alternate-hash-additive-averages--a). `-t` applies; `-o` does not. Mutually exclusive with `-n`. |
| `-d PATH` | With `-a`: outputs every computed key to `PATH` (sorted, one per line, as `KEY  path`) instead of anywhere else. Never appends — overwrites every run. Requires `-a`. |
| `-o N` | Max nilsimsa passes (default 50), `-n` only. |
| `-t N` | Number of worker threads, per frame: hashing (`-n`), averaging (`-a`), and writing the plane files (`-e`). Default: number of CPUs. |
| `-i SRC` | **Ingest:** read a y4m stream from `SRC` (a file, or `-` for stdin) into the new clip directory `<clipdir>`, creating or checking `specs.txt` from its header. Announces the destination directory before writing. See [Ingesting a y4m stream](#ingesting-a-y4m-stream--i). |
| *(several dirs)* | **Pooled sort:** `y4mstore DIR DIR ...` sorts several episodes (all frame clips, or all y/u/v component directories) into one sort file. See [Sorting several directories](#sorting-several-directories-at-once-episodes). |
| `-e SRC` | **Elementary planes:** read a y4m stream from `SRC` (a file, or `-` for stdin), write each frame's Y, U and V planes as separate files, then sort each plane directory (with `-n` or `-a`; default: no sort) and append the three lists into one. A directory that already holds `y/`, `u/` and `v/` planes is sorted the same way without `-e`. Announces the actual target directory before writing, including the `-elementary` rename when it triggers. See [Splitting into planes](#splitting-into-planes--e). |
| `-c` | **Check:** look for `specs.txt` in `<setdir>`; create it if missing (asking for the values), list it if present. Nothing is scanned, sorted or served. See [Checking specs.txt](#checking-specstxt--c). |
| `-s` | **Serve:** stream `<clipdir>` as a yuv4mpeg (y4m) stream. See [Serving as y4m](#serving-as-y4m). `-o`, `-t` don't apply. |
| `-f PATH` | Writes the sort file. Combined with `-s`, it writes the y4m stream to `PATH` instead. Never appends — prompts before overwriting. |
| `-p` | With `-i` or `-e` reading stdin (`-`): show a running count of frames read, so a slow upstream doesn't look like a hang. See [Progress on a piped stream](#progress-on-a-piped-stream--p). `-q` wins. |
| `-q` | Quiet — suppress progress bars and status messages. |
| `-h` | Show help. Printed to stdout, so it can be piped (`y4mstore -h \| less`); usage errors go to stderr with exit status 1. |

Threading only affects hashing wall-clock time: each digest is written to a fixed array slot by index, so output is identical regardless of `-t`.

## Directory layout and specs.txt

A *set* directory contains `specs.txt` and one subdirectory per clip:

```
myset/
  specs.txt
  clip01/   frame files...
  clip02/   frame files...
```

**Which directory is which is decided by what it contains, not by where a `specs.txt` happens to be.** A directory that holds frame files directly is a *clip*; its `specs.txt` is the one in its **parent**, whichever mode you name it in (`-s`, `-c`, or the sort). A directory with no frame files of its own is a *set*, and its own `specs.txt` applies.

- **A copy inside a clip is ignored** with a warning that it can be deleted. It is never streamed as a frame.
- **If the parent has none** but the clip itself has one, that copy is used, with a note suggesting you move it up to the set directory.
- **If neither exists**, you are asked for the values and the file is created in the **parent**, never inside the clip.

`specs.txt` describes the raw YUV common to every clip in the set, in five lines. The values use yuv4mpeg terms, so each maps directly onto a y4m header tag:

```
resolution=640x480
sampling=420
framerate=24000:1001
bitdepth=8
aspect=1:1
```

| Key | Valid values | y4m header |
|---|---|---|
| `resolution` | `WxH`, each side 1–65535. Example: `resolution=640x480` | `W640 H480` |
| `sampling` | `420`, `422`, `444`. Example: `sampling=420` | `C420jpeg`, `C422`, `C444` |
| `framerate` | A fraction `N:D` only. Examples: `24000:1001`, `24:1`, `25:1`, `30000:1001`, `60:1`. Decimals (`29.97`) and bare integers (`24`) are rejected; write whole rates as `24:1`. A `/` is accepted in place of `:` (`24000/1001`). Stored reduced (`48000:2002` becomes `24000:1001`). | `F24000:1001` |
| `bitdepth` | `8` or `10`. Example: `bitdepth=10` | no suffix for 8; `p10` suffix for 10 (`C420p10`, `C422p10`, `C444p10`) |
| `aspect` | `1:1` only (pixel aspect ratio). Never asked for. | `A1:1` |

- **Plane order:** frame files must hold planes as y4m does: all of Y, then all of U, then all of V (I420-style, not YV12).
- **10-bit:** each sample is a 16-bit little-endian word, as in y4m's `p10` colorspaces (`yuv420p10le`).
- **Square pixels:** footage is expected to be resampled to 1:1 pixels before it is stored, so `aspect=1:1` is the only accepted value; any other (`10:11`, `0:0`) is refused with a message. The header always carries `A1:1`.
- **Fixed header fields:** `specs.txt` doesn't carry interlacing, so the stream is always written as progressive (`Ip`). 4:2:0 is always tagged `420jpeg` (the y4m default chroma siting).

If `specs.txt` is missing and stdin is a terminal, y4mstore asks for the first four values one at a time (never `aspect`). Each question shows the valid values, an example line, and the y4m header token it produces, and re-asks on bad input. You can type either the bare value (`640x480`) or the whole example line (`resolution=640x480`). The answers are written to `specs.txt`, together with `aspect=1:1`, and it says so. If stdin is not a terminal it exits with an error instead of guessing. If an existing, otherwise valid `specs.txt` has no `aspect` line, y4mstore appends `aspect=1:1` and prints `added aspect=1:1 to <path>` (shown even with `-q`, since it edits your file); if the file can't be written it warns and carries on with 1:1. An existing `specs.txt` that is malformed is reported and left untouched: this includes files from earlier versions that used `i420`, `yv12`, `gray`, `bitdepth=12` and so on, or that lack `bitdepth`. Blank lines and `#` comments are ignored, and unknown keys are ignored with a warning.

### Checking specs.txt (`-c`)

`y4mstore -c <dir>` only deals with `specs.txt`; it doesn't scan, sort or serve anything. Name the set directory, or a clip directory (in which case it checks the set directory above it and says so).

- **Not found:** it asks for the four values (with the same prompts and examples as above), writes `specs.txt` with `aspect=1:1`, says `created <path>`, and then lists it.
- **Found:** it lists the specs, each with the y4m header token it produces, followed by the frame size and the complete header line:

```
myset/specs.txt
  resolution  640x480        -> W640 H480
  sampling    420            -> C420jpeg     (4:2:0 - chroma half width, half height)
  framerate   24000:1001     -> F24000:1001  (23.976 fps)
  bitdepth    8              -> 1 byte per sample (no suffix)
  aspect      1:1            -> A1:1
  frame size  460800 bytes
  y4m header  YUV4MPEG2 W640 H480 F24000:1001 Ip A1:1 C420jpeg
```

- The listing goes to stdout and is shown even with `-q`; status messages go to stderr.
- An existing `specs.txt` with no `aspect` line gets `aspect=1:1` added first, as usual. An invalid one is reported and left untouched.
- Exit status is 0 when `specs.txt` exists and is valid, 1 otherwise; it can't be combined with `-s`.

`specs.txt` itself is never treated as a frame. On a successful load y4mstore prints one line (suppressed by `-q`) naming the file it read, e.g. `specs (myset/specs.txt): 640x480 420 8-bit, 24000:1001 fps, 460800 bytes/frame`. The bytes-per-frame figure is worth checking against your actual file sizes.

Reading a file list from stdin (`-`) or passing a single file involves no directory layout, so no `specs.txt` is used there.

**Current status:** `-s` (below) uses the specs to build the y4m header. The similarity sort does not use them yet. Naming a clip directory sorts just that clip (one episode, one batch); naming a set directory still weighs all its clips together in one batch, as before.

## Ingesting a y4m stream (`-i`)

`-i` is the reverse of `-s`: it reads a yuv4mpeg stream and turns it into a clip directory of raw frame files, with `specs.txt` filled in from the stream's header. The source is the option's argument (`-` is stdin), and the positional argument is the **new clip directory**:

```
ffmpeg -i ep1.mkv -pix_fmt yuv420p -f yuv4mpegpipe - | y4mstore -i - myset/ep1
y4mstore -i clip.y4m myset/clip01
```

Each frame becomes one raw file, `yiv-000001.yuv`, `yiv-000002.yuv`, and so on, in the order they arrive. The bytes are copied through unchanged, since a y4m frame already is raw planar Y, U, V in the layout the other modes read. The result works straight away with `-s`, `-a`, `-c` and the sort; ingesting and then serving gives back pixel-identical video.

**`specs.txt` is checked before anything is written.** The header gives the resolution, sampling, frame rate and bit depth, and the clip's parent (the set directory) is where `specs.txt` lives:

- **Missing:** it is created from the header, with `aspect=1:1`.
- **Present and the same:** y4mstore says it matches and carries on. "Same" compares the values, so `48000:2002` equals `24000:1001`, and a slash-form frame rate is fine. If the file lacks the `aspect` line it is appended, as elsewhere.
- **Present and different:** it stops, lists each differing field (`specs.txt 64x48, stream 32x32`), and writes nothing: no frames, no clip directory, no changes to `specs.txt`.
- **Present but malformed:** it stops and leaves the file alone.

**What it accepts.** The colorspaces the tool can store: `420` (also `420jpeg`, `420mpeg2`, `420paldv`, and no `C` tag at all), `422`, `444`, and the 10-bit forms `420p10`, `422p10`, `444p10`. It refuses, before writing anything, anything else (`mono`, `411`, `444alpha`, 9/12/14/16-bit), a stream with no usable frame rate, and non-square pixels. This matches the rule that footage is resampled to 1:1 before it is stored. Unknown aspect (`A0:0`) is accepted and treated as square. Extra header tags such as ffmpeg's `XYSCSS` are ignored. Some things are noted rather than refused: the 420 chroma siting variants are stored as plain `420`, and an interlaced flag is not recorded, because `specs.txt` doesn't carry either. ffmpeg may need `-strict -1` to write 10-bit y4m (it did for 4:2:2 here).

- **New or empty destination only.** If the clip directory already holds anything, or is a file, it refuses before touching `specs.txt`. Frames are never mixed into an existing clip.
- **Damaged streams.** If the stream ends inside a frame, the partial frame is removed, the complete frames before it are kept, and the exit status is 1. A stream with a header but no frames, or a clip directory this run created but never filled, leaves nothing behind. A bogus header claiming an enormous frame size fails cleanly rather than allocating it.
- **Safe to run in parallel.** Several ingests into one set at the same time are fine: `specs.txt` is created atomically, so exactly one run writes it and the others see it and compare. When they disagree, whichever stream wrote `specs.txt` first is the one that proceeds.
- **Safe messages.** Text quoted from the stream in error messages is sanitized, so a damaged stream can't put control characters or escape sequences on your terminal.
- **Terminal guard.** `-i -` refuses to read from a terminal instead of waiting on it.
- **Progress:** a progress bar is shown for a file source (its size is known), and not for a pipe. `-q` silences the status lines.

Not combinable with `-a`, `-s` or `-c`. Measured on one core with a warm cache, 600 frames of 720x480 (about 311 MB) ingest in roughly 0.2 s, and the frames' checksum matched ffmpeg's own raw decode exactly.

### Progress on a piped stream (`-p`)

A y4m file has a known size, so `-i` and `-e` show a progress bar for one. A **piped** stream has no known size, so there is no bar, and a long `ffmpeg ... | y4mstore -i -` would otherwise print nothing until it finished. `-p` adds a running count of the frames read:

```
ffmpeg -i ep1.mkv -pix_fmt yuv420p -f yuv4mpegpipe - | y4mstore -p -i - myset/ep1
y4mstore -s myset/ep1 | y4mstore -p -e - myset/ep1
```

```
y4mstore: ingesting: 1234 frames, 532.1 MB, 41.2 frames/s, 0:30 elapsed
```

- **On a terminal** it is one line, redrawn about four times a second. **If stderr is redirected** (a log file) it is a plain line every couple of seconds instead, so a log doesn't fill with carriage-return redraws.
- **It appears immediately** with `0 frames`, before the first frame arrives, so you can tell it is waiting on upstream rather than hung. If the count stops rising, upstream has stalled.
- **Errors never get glued onto it:** any message starts on a fresh line.
- **It is opt-in** because upstream tools such as ffmpeg draw their own `\r` status line on the same terminal, and two programs redrawing one line garble each other. (Silence ffmpeg's with `-v error` or `-nostats` if you want both.)
- **ffmpeg's own logging goes to stderr by default, not the pipe,** so it doesn't corrupt the y4m stream on stdout either way. Adding `-loglevel quiet` (or `-hide_banner -v error`) to the ffmpeg command is still worth doing, though: it keeps ffmpeg's banner and per-frame stats out of your terminal entirely, so the only thing you see there is y4mstore's own progress line, with nothing to visually untangle from it.
- **`-q` wins** over `-p`. A **file** source keeps its progress bar and gets no counter. In every other mode `-p` is accepted and ignored. It only ever touches the reading thread, so it has no effect on the output or on threading.

## Splitting into planes (`-e`)

`-e` takes a y4m file or stream and splits every frame into its three components, then sorts each component's files and appends the results into one sort list. The input is a **y4m file or stream and nothing else**; to start from a clip directory of frames, pipe another instance in:

```
y4mstore -s myset/ep1 | y4mstore -e - myset/ep1
y4mstore -e clip.y4m myset/ep1
ffmpeg -i ep1.mkv -pix_fmt yuv420p -f yuv4mpegpipe - | y4mstore -e - myset/ep1 > all-planes.txt
```

A directory given as the source is refused with exactly that hint.

**Threading is per frame, as in `-a`.** A y4m stream can only be read in order, so the calling thread reads one whole frame at a time and hands it to one of the `-t` workers, which writes that frame's `y`, `u` and `v` files; a worker takes the next frame as soon as it finishes its last. The sort phase uses the same per-file workers as the other modes. Output is identical for any `-t`. A frame is only handed over once it has been read completely, so a truncated stream can never leave a frame with only some of its planes. Memory is bounded: a small fixed set of frame buffers is recycled (at most `-t` + 2 frames, and about 256 MiB of them), so the reader waits when the workers fall behind. The buffers grow only as data actually arrives, so a damaged header claiming a huge frame size can't cause a huge allocation.

**The planes.** For a subdirectory named `myset/ep1`:

```
myset/ep1/y/yuv-0000001.y     myset/ep1/u/yuv-0000001.u     myset/ep1/v/yuv-0000001.v
myset/ep1/y/yuv-0000002.y     ...
```

Each file is one plane, copied through unchanged: the Y plane is width x height samples, and the U and V planes have the size the sampling gives (86400 bytes for 4:2:0 at 720x480; twice that at 10-bit).

**The `-elementary` rule.** If the named subdirectory already exists and holds a `.yuv` file (a whole-frame clip, for instance from `-i`, which is exactly what you get when you pipe `-s myset/ep1` into `-e ... myset/ep1`), the planes go to `myset/ep1-elementary` instead, and y4mstore says so (even with `-q`). The original clip is untouched. In all cases the `y`, `u` and `v` directories must be new or empty, otherwise it refuses before writing anything; it never creates `-elementary-elementary`.

**`specs.txt`** is handled exactly as in `-i`: it lives in the set directory (the parent of the subdirectory) and is created from the stream header if missing, must match if present, and otherwise nothing is written.

**The sort.** Once the planes are written, the selected sorter runs on the `y` directory, then `u`, then `v`, and the three sorted lists are appended into one list, written to stdout or to `-f FILE`:

- **Default: no sort.** Files are listed in natural filename order, per plane; no hashing, no duplicate detection.
- **With `-n`: nilsimsa.** Each plane directory gets the normal greedy + 2-opt similarity sort. `-o` and `-t` apply to each plane. Files are listed in frame order before sorting, so the run is deterministic.
- **With `-a`: additive averages.** Each plane file is a single plane, so its key is that plane's own rounded average as **4 digits** (0 to 255 at 8-bit, up to 1023 at 10-bit), sorted ascending, with ties in filename order. Summary lines are labelled per plane (`[y] 41 distinct keys, ...`). A single plane average is a coarse key: expect many frames to share one, especially in the U and V planes, where averages cluster near mid-grey. The per-plane summary shows the distinct-key count straight away.

`-e` can be combined with `-n` or `-a` (not both) but not with `-i`, `-s` or `-c`.

**Sorting planes that already exist (no `-e` needed).** If the directory you name already holds `y/`, `u/` and `v/` planes, which is exactly what `-e` leaves behind, y4mstore recognises it on its own. Run it without `-e` and it does the second half of the job: the selected sorter runs on each component directory and the three lists are appended into one, y then u then v.

```
y4mstore myset/ep1-elementary                  # nilsimsa on each plane directory
y4mstore -a -f all-planes.txt myset/ep1-elementary    # additive averages, one combined file
```

The list is identical to the one `-e` printed when it created the planes (checked for both sorters, all formats, and at 600-frame scale), so you can re-sort later, or switch sorters, without splitting again.

- **How it is recognised:** `y/`, `u/` and `v/` exist and hold `.y`, `.u` and `.v` files respectively (case-insensitive), and the directory has no whole-frame `.yuv` files of its own. A directory with `.yuv` files is a clip whatever else it contains, and a set that merely has clips *named* `y`, `u` and `v` is not mistaken for one. Anything that doesn't match (a missing `v/`, no `.v` files) behaves exactly as it did before this feature.
- **Specs** come from the parent set directory, like a clip's. The additive sorter needs them for the plane sizes, so with `-a` a missing `specs.txt` is asked for (and created in the parent, never inside the component directory) or refused when there is no terminal. Nilsimsa needs none, so it only reads one that is already there and never prompts. A `specs.txt` inside the component directory is ignored with a warning when the parent has one, and used, with a note, when it doesn't.
- **Files** are read in natural filename order, so the result doesn't depend on directory order. If the three directories hold different numbers of files it warns and sorts each as it is.
- **Nothing is written unless all three sorts succeed.** The lists are collected in memory first, so with `-a` a plane of the wrong size (say, specs that don't match) exits with status 1 and leaves an existing `-f` file untouched.
- **Announced:** it says what it detected (silenced by `-q`). `-t`, `-o`, `-q` and `-f` apply as in `-e`.
- **The other modes agree:** `-e` into a directory that already holds all three planes refuses and tells you to run without `-e`; `-s` serves one back as a y4m stream (below); `-c` checks the parent's `specs.txt`.

**If the stream breaks.** If it ends inside a frame, that frame is never written at all, the complete frames before it are kept, no sort list is written, and the exit status is 1. If a write fails (a full disk, say), the failing frame's files are removed, the run stops without hanging, and every frame on disk has all three of its planes. A stream with no frames leaves no clip or plane directories behind (only the `specs.txt` made from its valid header, as with `-i`). The same header checks as `-i` apply, so unsupported colorspaces and non-square pixels are refused before anything is written, and text quoted from the stream in messages is sanitized. A `-f` file is only opened after the split succeeds, so an existing one survives a failed run.

Checked by comparing every plane file byte-for-byte with planes extracted independently from ffmpeg's own raw output (4:2:0, 4:2:2 10-bit, 4:4:4 and odd dimensions), by comparing each block of the combined list against the sorter run directly on that plane (and against a numpy reference for the additive mode), and by reassembling a 600-frame split into checksums identical to the original frames.

## Sorting several directories at once (episodes)

For sorting only, `y4mstore` accepts more than one directory and produces **one sort file** for all of them. The usual case is several episodes:

```
y4mstore -f all-episodes.txt set/ep1 set/ep2 set/ep3          # frame clips, nilsimsa
y4mstore -a -f all-episodes.txt set/ep*-elementary            # y/u/v component directories, additive
```

Directories are taken **in the order given**, so a shell glob gives episode order.

**One pooled sort, not one per episode.** All the directories' files are gathered and sorted together, so frames that recur across episodes (an opening sequence, a repeated shot) are grouped next to each other. With two copies of an episode, for example, the same frame from each copy ends up adjacent. Exact duplicates are counted over the whole pool too. The list is identical to what the plain sorter gives on the same combined list of files.

- **Frame clips:** every frame file from every directory is pooled and sorted once into one list.
- **Component directories:** all the `y` planes are pooled and sorted, then all the `u` planes, then all the `v` planes, and the three lists are appended: all the y's, then the u's, then the v's, as one sort file.

**Sanity checks, done before anything is written** (a failed check writes nothing and leaves an existing `-f` file alone):

- **Same type.** Every directory must be a whole-frame clip, or every directory must be a y/u/v component directory. Mixed lists are refused, with each directory's type shown. That includes a clip together with its own `-elementary` planes.
- **A real clip or component directory.** A set directory (clips inside, no frames of its own), an empty directory, or a half-made component directory (say, no `v/`) is refused. Name the clip directories instead.
- **Each directory once.** Naming a directory twice, in different spellings or through a symlink, is refused, since its files would be listed twice.
- **Existing directories,** and no `-` (a stdin list) mixed in.
- **Same format, for `-a`.** Additive keys from different resolutions, samplings or bit depths aren't comparable, so `-a` requires them to match (frame *rate* may differ). Nilsimsa works on any bytes, so it only warns.
- **Specs** are read from each directory's set, once per set. A component directory with different numbers of `y`, `u` and `v` files is warned about by name, and still sorted.

`-t`, `-o`, `-q` and `-f` apply as usual. nilsimsa's ordering cost grows roughly with the *square* of the pooled file count, so three episodes cost about nine times one; on 1,200 frames of 720x480 it took about 4 s against 0.25 s for `-a`.

Other modes take exactly one directory. `-s`, `-e`, `-i` and `-c` given several now say so and stop, instead of quietly ignoring the rest. One directory behaves exactly as before.

## Alternate hash: additive averages (`-a`)

An experimental alternative to the nilsimsa similarity sort. It needs no similarity hash and none of the greedy / 2-opt passes: each frame is reduced to one number, and the frames are sorted by it once.

1. **Split** each frame into its Y, U and V planes, using the resolution, sampling and bit depth from `specs.txt`.
2. **Add** every sample in a plane, then **divide** by that plane's own sample count. For a 720x480 frame the Y plane has 345600 samples; each 4:2:0 chroma plane has 86400. Averages are rounded to the nearest whole number (half rounds up).
3. **Format** each average as 4 digits (room for 10-bit values up to 1023) and join them, Y first: Y=250, U=127, V=36 gives `025001270036`.
4. **Sort** the frames numerically by that 12-digit key, ascending. Y is the most significant part, so frames are ordered by average luma first, then U, then V.

```
y4mstore -a -f list.txt myset/clip01
tar -cf - --no-recursion -T list.txt | zstd -T2 -13 --long > clip01.tar.zst
```

- **Ties are not duplicates.** Moving a pixel doesn't change any average, so different pictures can share a key. Frames with equal keys are ordered by filename (natural order), which keeps runs deterministic. Exact duplicates are still a separate thing. The run summary reports how many distinct keys there were and how many frames share a key; a lot of sharing means the key doesn't discriminate much on your footage.
- **This method cannot detect exact duplicates, by its nature.** Averaging is lossy and many-to-one: it collapses every frame down to three numbers, and two frames with the same average are not necessarily identical, while two frames that differ only slightly (or are byte-for-byte the same) will land on the same key or on adjacent ones, indistinguishably. A tie here is never proof of duplication, and the absence of a tie is never proof that two frames differ. If exact-duplicate detection is what you actually need, that has to come from something else (a real hash of the frame's bytes, or zpaq's own fragment-level dedup at the compression stage) — this mode isn't it, and isn't trying to be.
- **One frame per file.** Every file must be exactly one frame at the specs' frame size. If any isn't (a stray `.DS_Store`, a wrong resolution in `specs.txt`), y4mstore lists the offenders and exits without writing anything.
- **Threading is per frame.** Every frame is independent, so each worker in the `-t` pool claims the next single frame as soon as it finishes its last. One OS thread per frame would cost more than the work; this gets the same per-frame parallelism with balanced load. Output is identical for any `-t`.
- **Same layout rules as the other modes.** Name a clip directory and the specs come from its parent. Naming a set directory keys every clip in it together into one list. It needs a directory, not `-` or a single file.
- **`-d PATH` outputs every key to a log file,** sorted, one per line, as `KEY  path`, instead of anywhere else, so a long list doesn't scroll past the progress display. Never appends — each run overwrites `PATH` fully, and if it already exists you're asked before it's replaced. Requires `-a`.
- **Not a similarity measure.** Two frames with the same three averages can look nothing alike, and two very similar frames can straddle a rounding boundary. Whether it helps compression is an empirical question; compare against directory order and a `shuf`ed list.

**Measured (one core, page cache warm):** 400 frames of 720x480 4:2:0 took 0.11 s with `-a` (about 1.9 GB/s) versus 1.87 s for the nilsimsa mode on the same frames. Every thread count from 1 to 64 produced identical output. Correctness was checked against an independent numpy implementation across 4:2:0 / 4:2:2 / 4:4:4, 8 and 10-bit, odd dimensions, rounding edges and deliberate ties, and ThreadSanitizer found no races.

## Serving as y4m

`-s` turns a clip directory of raw frame files into a [yuv4mpeg](https://linux.die.net/man/5/yuv4mpeg) stream: one header line built from `specs.txt`, then `FRAME` plus the raw bytes for each frame. Name the **clip** directory; its `specs.txt` is read from the parent set directory, as described above. If it isn't there yet you are asked for it, and the file is created in the parent.

```
y4mstore -s myset/clip01 | ffplay -                        # watch it
y4mstore -s myset/clip01 | mpv -
y4mstore -s myset/clip01 | ffmpeg -i - -c:v ffv1 clip01.mkv   # convert
y4mstore -s myset/clip01 > clip01.y4m                      # save
y4mstore -s -f clip01.y4m myset/clip01                     # save, with a progress bar
```

For a 24000:1001, 640x480, 4:2:0, 8-bit set the stream starts with:

```
YUV4MPEG2 W640 H480 F24000:1001 Ip A1:1 C420jpeg
```

- **Frame order** is filename order, with digit runs compared as numbers (`f2` before `f10`). Dotfiles are skipped, and subdirectories inside the clip directory are ignored with a warning.
- **One frame per file, or several.** Each file must be a whole number of frames at the specs' frame size; a file holding three frames is served as three. The summary line shows `N frames from M files`, so a wrong resolution that happens to divide evenly stands out.
- **Nothing is written on a mismatch.** Every file is checked first; if any size doesn't match, y4mstore lists the offenders and exits without producing output.
- **Safeguards:** it refuses to write a video stream to a terminal, and names the fix if you point it at a set directory (one with no frame files of its own) instead of a clip.
- **Closing the player** (or `head`) stops it quietly with exit status 0.

Since the raw frames are already in y4m plane order, the bytes are copied through untouched.


### Serving a 3-component directory back to y4m

`-s` also takes a directory of `y/`, `u/` and `v/` planes (what `-e` writes) and streams it back out as a y4m stream. This completes the round trip: y4m in, planes on disk, y4m out.

```
y4mstore -s myset/ep1-elementary | ffplay -
y4mstore -s -f ep1.y4m myset/ep1-elementary
y4mstore -s myset/ep1-elementary | y4mstore -i - myset/ep1-restored     # back to whole frames
```

Each frame is its `y`, `u` and `v` file one after another, under a header built from `specs.txt` (in the parent set directory, as for a clip). The stream is **byte-identical to serving the same clip as whole frames** with `-s`, and decodes in ffmpeg to exactly the original video (checked for 4:2:0, 4:2:2 10-bit, 4:4:4 and odd dimensions, and for a 600-frame 720x480 run).

- **The three planes of a frame are paired by filename stem**, not by position: `yuv-0000007.y`, `.u` and `.v` make one frame. So a missing or stray file can never quietly mix planes from different frames. If any frame lacks a plane, y4mstore names it (`frame 'yuv-0000012': found y v, missing u`) and writes nothing.
- **Every plane file must be exactly the plane size** the specs give, otherwise it lists the offenders and writes nothing.
- **Order** is natural filename order by stem (`f2` before `f10`), the same as whole-frame `-s`. The stems don't have to be the `yuv-0000001` names `-e` writes, as long as `y/`, `u/` and `v/` agree.
- The recognition rule, the safeguards (no terminal, clean stop when the reader quits, `-f FILE`, `-q`) and the specs lookup are the same as described above. A directory that also holds `.yuv` files is a clip and is served as whole frames.
- **The header is regenerated from `specs.txt`,** so extra header tags in an original file (ffmpeg's `XYSCSS`, `XCOLORRANGE`) are not carried through; the video data is. This is the same for every `-s` and `-i`.

## Intended pairing: zpaq + zpaqfuse, and a caveat for HDDs

The point of pooling several episodes or a similar-content video set into
one y4mstore store, rather than keeping each as its own archive, is
deduplication: pool them, sort the pooled frames, and hand that order to
zpaq. Similar or repeated content across episodes, recurring backgrounds,
openings, reused cels, only has a chance of landing in the same zpaq block
if the sort actually puts it next to itself; sorted within one episode at a
time, it never gets that chance. This is the motivation for y4mstore's
multi-directory pooled sort (see "Sorting several directories at once").

y4mstore is meant to be paired with [zpaq](https://mattmahoney.net/dc/zpaq.html)
(patched to add `-order`, so it takes a y4mstore sort list directly, without
a tar) and zpaqfuse to mount the resulting archive read-only. y4mstore stores and sorts the frames; zpaq compresses them into one
small archive; zpaqfuse mounts that archive so `-s` (or anything else) can
read individual clips or planes straight out of it, without extracting first.
That combination is what it's designed for: serving video off an SSD-backed
store, where the archive can be mounted and read from directly.

**On a spinning disk (HDD), don't mount and serve from the archive.** Reading
a clip back out means jumping to wherever its blocks live inside the zpaq
archive, and repeating that for every clip or plane you serve turns into a
lot of random seeks, exactly what an HDD is worst at. On an SSD there's no
seek penalty, so this isn't a concern there.

For HDD storage, use this pipeline for archiving only, not for serving:

```
y4mstore myset/ep1 | tar -T- --no-recursion -cf - | xz -T1 --lzma2=preset=6,dict=32MiB > ep1.tar.xz
```

A single sequential `.tar.xz` stream reads back with one long sequential pass
instead of scattered seeks, which suits an HDD far better than mounting a
block-structured archive on it. Extract it back to disk fully before serving
from it, rather than mounting and reading from the compressed stream directly.

## Duplicate detection

Alongside `-n`'s sort, y4mstore tracks exact duplicates — files byte-identical over their full contents — using a 64-bit FNV-1a fingerprint plus a length check computed in the same read pass. This is a byproduct of `-n`'s hashing specifically: it does not run by default, and it does not run with `-a`. At the end of a run it prints a one-line summary to stderr (suppressed by `-q`):

```
y4mstore: 4 duplicate files found in 2 sets of identical content
```

This is a count only; nothing extra goes to stdout, so piping the sort order is unaffected.

This is a separate mechanism from the `-a` additive-average sort. This duplicate count comes from a fingerprint over each file's actual bytes, so it can tell true duplicates apart from lookalikes. The additive average, by contrast, reduces a frame to three lossy numbers and cannot make that distinction: two frames sharing a key are not necessarily duplicates, and two duplicates are not guaranteed to land on the same key. Nothing from `-a` feeds into, or is validated by, this count.

## Notes for future work

- **Single-episode results vs. FFV1:** raw YUV, FFV1, and zpaq (via zpaq-order) sizes for one episode, to post here once run.
- **3-part pooled episode results vs. FFV1:** the same comparison, pooling three episodes into one sort and archive, to see how much cross-episode redundancy zpaq-order can actually capture.

## License

MIT. See `LICENSE` and the header comment in `y4mstore.c`. y4mstore derives from nilsort, so both copyright lines are retained.

The nilsimsa algorithm itself (the TRAN/POPC-derived tables, `tran3` mixing function, and digest construction) is a faithful port of the public reference implementation originally written by cmeclax, based on Damiani et al. 2004, "An Open Digest-based Technique for Spam Detection." Threading, greedy + 2-opt ordering, progress reporting, and the POPCNT-based comparison are inherited from nilsort.
