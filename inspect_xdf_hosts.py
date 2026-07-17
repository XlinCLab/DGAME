import sys

from pyxdf import load_xdf

from utils.xdf_utils import get_stream_hostname

if len(sys.argv) != 2:
    print("Usage: python inspect_xdf_hosts.py /path/to/file.xdf")
    sys.exit(1)

xdf_file = sys.argv[1]
streams, _header = load_xdf(xdf_file, verbose=False)

for stream in streams:
    info = stream["info"]
    name = info.get("name", [""])[0]
    stype = info.get("type", [""])[0]
    hostname = get_stream_hostname(stream)
    n_channels = info.get("channel_count", [""])[0]
    srate = info.get("nominal_srate", [""])[0]
    stream_id = info.get("stream_id", "?")
    segments = info.get("segments", [])
    clock_segments = info.get("clock_segments", [])
    print(f"name={name!r:25} type={stype!r:12} hostname={hostname!r:15} channels={n_channels} srate={srate} stream_id={stream_id}")
    print(f"    segments ({len(segments)}): {segments}")
    print(f"    clock_segments ({len(clock_segments)}): {clock_segments}")
    if len(segments) > 1:
        print(f"    -> REAL data gap(s): {len(segments) - 1} gap(s) between recorded samples")
    elif segments != clock_segments:
        print("    -> mismatch is only in clock-reset detection, no gap in recorded samples")
