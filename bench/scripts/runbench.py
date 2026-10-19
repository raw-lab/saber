#!/usr/bin/env python3
"""Run a command N times and record wall time, CPU time and peak RSS.

Uses os.wait4() so the numbers come from the kernel's rusage for the child
(and its own children are NOT included -- every tool here is single-process).

Usage:
  runbench.py --label NAME --reps 3 --out results.jsonl [--stdout FILE] -- cmd args...

Each repetition appends one JSON object to --out:
  {label, rep, cmd, wall_s, user_s, sys_s, cpu_s, maxrss_kb, exit, stdout_bytes,
   timestamp, cpu}
`cpu` (model name and last-level cache size from /proc/cpuinfo) identifies the
machine, so tables can be checked to compare runs from one host only.
The command's stdout is written to --stdout (default /dev/null) and stderr to
<out>.stderr/<label>.<rep>.log so failures are inspectable.
"""
import argparse, json, os, shlex, subprocess, sys, time

def cpu_id():
    model = cache = "?"
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name") and model == "?":
                    model = line.split(":", 1)[1].strip()
                elif line.startswith("cache size") and cache == "?":
                    cache = line.split(":", 1)[1].strip()
    except OSError:
        pass
    return f"{model} (cache {cache})"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stdout", default=os.devnull)
    ap.add_argument("--drop-caches", action="store_true",
                    help="attempt 'echo 3 > /proc/sys/vm/drop_caches' before each rep (cold cache)")
    ap.add_argument("--env", action="append", default=[], help="KEY=VALUE for the child")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        sys.exit("no command")
    errdir = a.out + ".stderr"
    os.makedirs(errdir, exist_ok=True)
    env = dict(os.environ)
    for kv in a.env:
        k, v = kv.split("=", 1)
        env[k] = v
    for rep in range(a.reps):
        cold = False
        if a.drop_caches:
            try:
                subprocess.run(["sync"], check=False)
                with open("/proc/sys/vm/drop_caches", "w") as f:
                    f.write("3\n")
                cold = True
            except Exception as e:  # not permitted in most containers
                cold = False
        errpath = os.path.join(errdir, f"{a.label}.{rep}.log")
        with open(a.stdout, "wb") as so, open(errpath, "wb") as se:
            t0 = time.perf_counter()
            p = subprocess.Popen(cmd, stdout=so, stderr=se, env=env)
            _, status, ru = os.wait4(p.pid, 0)
            wall = time.perf_counter() - t0
        code = os.waitstatus_to_exitcode(status)
        try:
            sz = os.path.getsize(a.stdout) if a.stdout != os.devnull else 0
        except OSError:
            sz = 0
        rec = dict(label=a.label, rep=rep, cmd=" ".join(shlex.quote(c) for c in cmd),
                   wall_s=round(wall, 4), user_s=round(ru.ru_utime, 4), sys_s=round(ru.ru_stime, 4),
                   cpu_s=round(ru.ru_utime + ru.ru_stime, 4), maxrss_kb=ru.ru_maxrss,
                   exit=code, stdout_bytes=sz, cold_cache=cold,
                   timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"), cpu=cpu_id())
        with open(a.out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(json.dumps(rec), file=sys.stderr)
        if code != 0:
            sys.exit(code)

if __name__ == "__main__":
    main()
