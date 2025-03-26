#!/usr/bin/env python3

import subprocess
import time
import sys
import os
import platform
import argparse
import multiprocessing

def parse_circuit_output(output):
    """Parse the program output to extract the optimal circuit and its cost."""
    # Look for the circuit pattern in the output
    # The circuit is printed as pairs of lines: gates and cost
    circuits = []
    current_circuit = []
    current_cost = None
    
    for line in output.split('\n'):
        line = line.strip()
        # Skip empty lines, progress indicators, timing information, and log messages
        if (not line or 
            line.startswith('|') or 
            line.startswith('Looking for') or 
            line.startswith('Time:') or
            line.startswith('Searching for') or
            line.startswith('Iteration set size:') or
            line.startswith('Found')):
            continue
            
        # If line starts with "Cost"
        if line.startswith('Cost'):
            try:
                current_cost = int(line.split()[1])
                if current_circuit:
                    circuits.append(('\n'.join(current_circuit), current_cost))
                    current_circuit = []
                    current_cost = None
            except (IndexError, ValueError):
                continue
        # If line contains gates
        elif any(gate in line for gate in ['I', 'H', 'X', 'Y', 'Z', 'S', 'T', 'C']):
            current_circuit.append(line)
    
    if not circuits:
        print("\nDebug - No circuits found in output. Raw output:")
        print(output)
        return None, None
        
    # Find the circuit with lowest cost
    optimal_circuit, optimal_cost = min(circuits, key=lambda x: x[1])
    return optimal_circuit, optimal_cost

def compare_circuits(circuit1, circuit2):
    """Compare two circuits and return detailed differences."""
    if circuit1 == circuit2:
        return None
    
    # Split circuits into lines for comparison
    lines1 = circuit1.split('\n')
    lines2 = circuit2.split('\n')
    
    # Find the first line where they differ
    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        if line1 != line2:
            return f"First difference at line {i+1}:\nOld: {line1}\nNew: {line2}"
    
    # If one circuit is longer than the other
    if len(lines1) != len(lines2):
        return f"Different number of lines:\nOld: {len(lines1)} lines\nNew: {len(lines2)} lines"
    
    return "Unknown difference"

def run_program(program_path, args, num_times=1):
    """Run a program with given arguments and return its execution time and circuit output."""
    try:
        # First run to capture output
        result = subprocess.run([program_path] + args, 
                              capture_output=True, 
                              text=True)
        
        # Check if program executed successfully
        if result.returncode != 0:
            print(f"Error running {program_path}:")
            print(result.stderr)
            return None, None, None
            
        # Parse the circuit from the output
        circuit, cost = parse_circuit_output(result.stdout)
        if circuit is None:
            print(f"\nDebug - Failed to parse output from {program_path}")
            return None, None, None

        # Now run multiple times to get average runtime
        total_time = 0
        for _ in range(num_times):
            start_time = time.time()
            subprocess.run([program_path] + args, 
                         capture_output=True, 
                         text=True)
            total_time += time.time() - start_time

        avg_time = total_time / num_times
        return avg_time, circuit, cost
    except Exception as e:
        print(f"Error running {program_path}: {str(e)}")
        return None, None, None

def main():
    parser = argparse.ArgumentParser(description='Compare runtime of two MITMS program versions')
    parser.add_argument('circuit', help='The quantum circuit to synthesize')
    parser.add_argument('--threads', type=int, 
                       default=multiprocessing.cpu_count(),
                       help='Number of threads to use (default: number of CPU cores)')
    parser.add_argument('--times', type=int,
                       default=1,
                       help='Number of times to run each program for averaging (default: 1)')
    args = parser.parse_args()

    # Get the platform-specific paths
    current_platform = platform.system().lower()
    old_program = f"old/{current_platform}/mitms"
    new_program = "build/mitms"

    # Check if both programs exist
    if not os.path.exists(old_program):
        print(f"Error: Old program not found at {old_program}")
        sys.exit(1)
    if not os.path.exists(new_program):
        print(f"Error: New program not found at {new_program}")
        sys.exit(1)

    # Prepare arguments
    program_args = [
        "-threads", str(args.threads),
        "-no-equiv-checks",
        "-early-stop",
        args.circuit
    ]

    print(f"Running old version ({old_program})...")
    print(f"Arguments: {' '.join(program_args)}")
    print(f"Running {args.times} times for averaging")
    old_time, old_circuit, old_cost = run_program(old_program, program_args, args.times)
    
    print(f"\nRunning new version ({new_program})...")
    print(f"Arguments: {' '.join(program_args)}")
    print(f"Running {args.times} times for averaging")
    new_time, new_circuit, new_cost = run_program(new_program, program_args, args.times)

    if old_time is not None and new_time is not None:
        print("\nResults:")
        print(f"Old version: {old_time:.3f} seconds (average of {args.times} runs)")
        print(f"New version: {new_time:.3f} seconds (average of {args.times} runs)")
        print(f"Speedup: {old_time/new_time:.2f}x")
        
        print("\nCircuit Comparison:")
        if old_circuit and new_circuit:
            print(f"Old version optimal circuit (cost: {old_cost}):")
            print(old_circuit)
            print(f"\nNew version optimal circuit (cost: {new_cost}):")
            print(new_circuit)
            
            if old_circuit == new_circuit:
                print("\n✓ Circuits are identical")
            else:
                print("\n✗ Circuits are different!")
                diff = compare_circuits(old_circuit, new_circuit)
                if diff:
                    print("\nDetailed difference:")
                    print(diff)
        else:
            print("Could not parse circuits from output")
    else:
        print("\nComparison failed due to errors. See above for details.")

if __name__ == "__main__":
    main() 