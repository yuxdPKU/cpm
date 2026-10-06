#!/usr/bin/env python3
"""Run a command in an owned process group; always clean and verify the group."""
import argparse
import os
import signal
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--timeout', type=float, required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or not 0 < args.timeout < float('inf'):
        parser.error('a command and finite positive timeout are required')
    process = subprocess.Popen(command, start_new_session=True)
    interrupted = False

    def interrupt(signum, _frame):
        nonlocal interrupted
        interrupted = True
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, interrupt)
    signal.signal(signal.SIGINT, interrupt)
    code = 1
    try:
        code = process.wait(timeout=args.timeout)
    except subprocess.TimeoutExpired:
        print('Guarded command timed out', file=sys.stderr)
        code = 124
    except KeyboardInterrupt:
        code = 130
    finally:
        # Cleanup itself cannot be interrupted halfway through.
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(process.pid, sig)
            except ProcessLookupError:
                break
            deadline = time.monotonic() + 3
            while time.monotonic() < deadline:
                process.poll()  # reap the direct child if needed
                try:
                    os.killpg(process.pid, 0)
                except ProcessLookupError:
                    break
                time.sleep(.05)
            else:
                continue
            break
        process.wait()
        try:
            os.killpg(process.pid, 0)
        except ProcessLookupError:
            pass
        else:
            raise RuntimeError(f'cleanup verification failed for owned process group {process.pid}')
    return 130 if interrupted else (code if code >= 0 else 128 - code)


if __name__ == '__main__':
    sys.exit(main())
