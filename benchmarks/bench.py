#!/usr/bin/env python3
"""Run a command, report wall time and peak RSS (via wait4 rusage, which
captures the child's true high-water mark including all its threads)."""
import os, sys, time, resource, subprocess, json

def run(label, cmd, outfile=None, reps=1):
    best = None
    for _ in range(reps):
        out_fd = open(outfile, "wb") if outfile else subprocess.DEVNULL
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=out_fd, stderr=subprocess.DEVNULL)
        pid, status, ru = os.wait4(p.pid, 0)
        wall = time.perf_counter() - t0
        if outfile: out_fd.close()
        rss_mb = ru.ru_maxrss / 1024.0          # ru_maxrss is KB on Linux
        cpu = ru.ru_utime + ru.ru_stime
        rec = dict(label=label, wall=round(wall, 3), rss_mb=round(rss_mb, 1),
                   cpu=round(cpu, 3), status=status)
        if best is None or wall < best["wall"]:
            best = rec
    return best

if __name__ == "__main__":
    label = sys.argv[1]
    outfile = sys.argv[2] if sys.argv[2] != "-" else None
    reps = int(sys.argv[3])
    cmd = sys.argv[4:]
    r = run(label, cmd, outfile, reps)
    print(json.dumps(r))
