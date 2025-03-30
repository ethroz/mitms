#!/usr/bin/env python3

import subprocess
import sys
import argparse
from typing import List

def run_comparison(circuit: str, threads: int, times: int) -> bool:
    """Run the compare_runtime script with the given parameters."""
    cmd = ["./compare_runtime.py", circuit, "--threads", str(threads), "--times", str(times)]
    
    print(f"\nRunning comparison with {threads} threads...")
    result = subprocess.run(cmd)
    return result.returncode == 0

def main():
    parser = argparse.ArgumentParser(description='Compare MITMS program versions with different thread counts')
    parser.add_argument('circuit', help='The quantum circuit to synthesize')
    parser.add_argument('threads', nargs='+', type=int, help='List of thread counts to test')
    parser.add_argument('--times', type=int, default=1,
                       help='Number of times to run each program for averaging (default: 1)')
    args = parser.parse_args()

    # Validate arguments
    if args.times < 1:
        print("Error: Number of times to run each program must be at least 1")
        sys.exit(1)
    
    if not all(t > 0 for t in args.threads):
        print("Error: All thread counts must be positive integers")
        sys.exit(1)

    # Run comparison for each thread count
    success = True
    for threads in args.threads:
        if not run_comparison(args.circuit, threads, args.times):
            success = False
            print(f"\nComparison failed with {threads} threads")
            break

    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main() 