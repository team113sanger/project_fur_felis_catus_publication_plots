#!/bin/bash
# Download the feline reference genome from Ensembl
# Usage: ./download_feline_reference_genome.sh [OUTPUT_FILE]
set -euo pipefail

# CONSTANTS
DEFAULT_OUTPUT_FILE="$(pwd)/Felis_catus.Felis_catus_9.0.dna.toplevel.fa"
DOWNLOAD_URL="http://ftp.ensembl.org/pub/release-104/fasta/felis_catus/dna/Felis_catus.Felis_catus_9.0.dna.toplevel.fa.gz"
CHECKSUMS_URL="http://ftp.ensembl.org/pub/release-104/fasta/felis_catus/dna/CHECKSUMS"
FILE_BASENAME=$(basename "$DOWNLOAD_URL")
REQUIRED_PROGRAMS=("curl" "gunzip" "du" "df" "awk" "sum")
MINIMUM_DISK_SPACE_MB=4000  # Approximate minimum disk space required in MB
CURL_TIMEOUT=1800           # Timeout in seconds (30 minutes)
CURL_RETRY=3                # Number of retries
CURL_RETRY_DELAY=15         # Delay between retries in seconds
TEMP_DIR=""                 # Will be set to a temporary directory

# GLOBAL VARIABLES
OUTPUT_FILE=""
FORCE_OVERWRITE=false

print_usage() {
    echo "Usage: $0 [OPTIONS] [OUTPUT_FILE]"
    echo "Download the Felis_catus_9.0 genome from Ensembl and decompress it to a FASTA file."
    echo
    echo "Options:"
    echo "  -f, --force   Overwrite output file if it exists"
    echo "  -h, --help    Display this help message and exit"
    echo
    echo "Arguments:"
    echo "  OUTPUT_FILE   The output file path (default: ${DEFAULT_OUTPUT_FILE})"
    echo
    echo "Required programs:"
    for program in "${REQUIRED_PROGRAMS[@]}"; do
        echo "  $program"
    done
}

print_info() {
    # Print information messages in green to standard error
    echo -e "\e[32m[INFO] $1\e[0m" >&2
}

print_warning() {
    # Print warning messages in yellow to standard error
    echo -e "\e[33m[WARNING] $1\e[0m" >&2
}

print_error() {
    # Print error messages in red to standard error
    echo -e "\e[31m[ERROR] $1\e[0m" >&2
}

cleanup() {
    # Clean up temporary files
    if [[ -n "$TEMP_DIR" ]] && [[ -d "$TEMP_DIR" ]]; then
        print_info "Cleaning up temporary files..."
        rm -rf "$TEMP_DIR"
    fi
}

handle_error() {
    print_error "$1"
    cleanup
    exit 1
}

# Set up trap for script termination
trap cleanup EXIT
trap 'handle_error "Script interrupted by user. Exiting..."' INT TERM

check_disk_space() {
    local dir="$1"
    local required_mb="$2"
    
    # Get available disk space in KB and convert to MB
    local available_kb=$(df -P "$dir" | awk 'NR==2 {print $4}')
    local available_mb=$((available_kb / 1024))
    
    if (( available_mb < required_mb )); then
        handle_error "Not enough disk space in $dir. Required: ${required_mb}MB, Available: ${available_mb}MB"
    fi
    
    print_info "Sufficient disk space available: ${available_mb}MB"
}

check_dependencies() {
    for program in "${REQUIRED_PROGRAMS[@]}"; do
        if ! command -v "$program" &> /dev/null; then
            handle_error "$program is not installed. Please install it and try again."
        fi
    done
    print_info "All required dependencies are installed."
}

parse_args() {
    # Parse command-line arguments & update global variables
    local output_file=""
    
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -f|--force)
                FORCE_OVERWRITE=true
                ;;
            -h|--help)
                print_usage
                exit 0
                ;;
            -*)
                print_error "Unknown option: $1"
                print_usage
                exit 1
                ;;
            *)
                if [[ -z "$output_file" ]]; then
                    output_file="$1"
                else
                    print_error "Unknown argument: $1"
                    print_usage
                    exit 1
                fi
                ;;
        esac
        shift
    done
    
    OUTPUT_FILE="${output_file:-$DEFAULT_OUTPUT_FILE}"
    
    # Check that the output file name doesn't have a .gz extension
    if [[ "$OUTPUT_FILE" == *.gz ]]; then
        handle_error "The output file name should not have a .gz extension."
    fi
}

create_temp_dir() {
    # Create a temporary directory
    TEMP_DIR=$(mktemp -d 2>/dev/null || mktemp -d -t 'feline_genome_XXXXXX')
    if [[ ! -d "$TEMP_DIR" ]]; then
        handle_error "Failed to create temporary directory."
    fi
    print_info "Created temporary directory: $TEMP_DIR"
}

download_checksums() {
    local checksums_file="$1"
    local checksums_url="$2"
    
    print_info "Downloading CHECKSUMS file from Ensembl..."
    
    if ! curl --connect-timeout 30 \
              --max-time 60 \
              --retry 3 \
              --silent \
              --show-error \
              -L -o "$checksums_file" "$checksums_url"; then
        handle_error "Failed to download the CHECKSUMS file."
    fi
    
    # Check if the file was actually downloaded and has content
    if [[ ! -s "$checksums_file" ]]; then
        handle_error "Downloaded CHECKSUMS file is empty or does not exist."
    fi
    
    print_info "CHECKSUMS file downloaded successfully."
}

verify_checksum() {
    local file="$1"
    local filename="$2"
    local checksums_file="$3"
    
    print_info "Verifying file integrity using checksum..."
    
    # Look for the file entry in the CHECKSUMS file
    local checksum_entry=$(grep -E "[0-9]+ [0-9]+ $filename\$" "$checksums_file")
    
    if [[ -z "$checksum_entry" ]]; then
        handle_error "Could not find entry for $filename in CHECKSUMS file."
    fi
    
    # Parse the expected checksum and file size from the CHECKSUMS file
    local expected_checksum=$(echo "$checksum_entry" | awk '{print $1}')
    local expected_size_kb=$(echo "$checksum_entry" | awk '{print $2}')
    
    # Check file size (converting actual bytes to KB for comparison)
    local actual_size_bytes=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file")
    local actual_size_kb=$((actual_size_bytes / 1024))
    
    # Allow a small margin of error (±1KB) in the size comparison due to potential rounding
    if [[ "$actual_size_kb" -lt $((expected_size_kb - 1)) ]] || [[ "$actual_size_kb" -gt $((expected_size_kb + 1)) ]]; then
        handle_error "File size mismatch for $filename. Expected: ${expected_size_kb} KB, Got: ${actual_size_kb} KB (${actual_size_bytes} bytes)."
    fi
    
    print_info "File size verified: ${actual_size_kb} KB (${actual_size_bytes} bytes)"
    
    # Calculate the Unix sum checksum
    local actual_checksum=$(sum "$file" | awk '{print $1}')
    
    # STRING comparison to check if the checksums match
    if [ "$actual_checksum" != "$expected_checksum" ]; then
        handle_error "Checksum mismatch for $filename. Expected: $expected_checksum, Got: $actual_checksum."
    fi
    
    print_info "Checksum verified: $actual_checksum"
    print_info "File integrity confirmed."
}

download_genome() {
    local temp_file="$1"
    local download_url="$2"
    local temp_checksums_file="$TEMP_DIR/CHECKSUMS"
    
    # First download the CHECKSUMS file
    download_checksums "$temp_checksums_file" "$CHECKSUMS_URL"
    
    print_info "Downloading Felis_catus_9.0 genome from Ensembl..."
    print_info "URL: $download_url"
    print_info "This may take a while depending on your internet connection..."
    
    # Download with progress bar, timeout, and retry
    if ! curl --connect-timeout 30 \
              --max-time "$CURL_TIMEOUT" \
              --retry "$CURL_RETRY" \
              --retry-delay "$CURL_RETRY_DELAY" \
              --retry-max-time $((CURL_TIMEOUT * 2)) \
              -# -L -o "$temp_file" "$download_url"; then
        handle_error "Failed to download the genome file."
    fi
    
    # Check if the file was actually downloaded and has content
    if [[ ! -s "$temp_file" ]]; then
        handle_error "Downloaded file is empty or does not exist."
    fi
    
    local file_size=$(du -h "$temp_file" | cut -f1)
    print_info "Download complete. File size: $file_size"
    
    # Verify checksum and file size
    verify_checksum "$temp_file" "$FILE_BASENAME" "$temp_checksums_file"
}

decompress_genome() {
    local compressed_file="$1"
    local output_file="$2"
    local output_dir=$(dirname "$output_file")
    
    print_info "Decompressing the downloaded genome file..."
    
    # Ensure the output directory exists
    if [[ ! -d "$output_dir" ]]; then
        print_info "Creating output directory: $output_dir"
        mkdir -p "$output_dir" || handle_error "Failed to create output directory: $output_dir"
    fi
    
    # Decompress the file directly to the output location
    if ! gunzip -c "$compressed_file" > "$output_file"; then
        handle_error "Failed to decompress the genome file."
    fi
    
    # Verify the output file exists and has content
    if [[ ! -s "$output_file" ]]; then
        handle_error "Decompressed file is empty or does not exist."
    fi
    
    local file_size=$(du -h "$output_file" | cut -f1)
    print_info "Decompression complete. File size: $file_size"
}

main() {
    # Check if output file already exists
    if [[ -f "$OUTPUT_FILE" ]]; then
        if $FORCE_OVERWRITE; then
            print_warning "Output file already exists. Will overwrite: $OUTPUT_FILE"
        else
            handle_error "Output file already exists: $OUTPUT_FILE. Use -f or --force to overwrite."
        fi
    fi

    # Check if required programs are installed
    check_dependencies
    
    # Create temporary directory
    create_temp_dir
    
    # Check disk space in output directory and temp directory
    local output_dir=$(dirname "$OUTPUT_FILE")
    check_disk_space "$output_dir" "$MINIMUM_DISK_SPACE_MB"
    check_disk_space "$TEMP_DIR" "$MINIMUM_DISK_SPACE_MB"
    
    # Define temporary file for compressed download
    local temp_compressed="$TEMP_DIR/$(basename "$DOWNLOAD_URL")"
    
    # Download the genome file
    download_genome "$temp_compressed" "$DOWNLOAD_URL"
    
    # Decompress the genome file
    decompress_genome "$temp_compressed" "$OUTPUT_FILE"
    
    print_info "Process completed successfully."
    print_info "Genome saved to: $OUTPUT_FILE"
}

# RUN THE SCRIPT
parse_args "$@"
main