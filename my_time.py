import argparse
import subprocess
import time

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Measure the execution time of a command.')
parser.add_argument('command', nargs=argparse.REMAINDER, help='The command to execute and time')
args = parser.parse_args()

# Ensure a command is provided
if not args.command:
    print("Error: No command provided.")
    exit(1)

# Start timing
start_time = time.time()

# Execute the command
try:
    subprocess.run(args.command)
except Exception as e:
    print(f"Error executing command: {e}")
    exit(1)

# End timing
end_time = time.time()

# Calculate elapsed time
elapsed_time = end_time - start_time

# Output the execution time
print(f"Elapsed time: {elapsed_time:.3f} seconds") 