#!/bin/bash
# Test script for --connected_blocks feature
# Tests kaffpa with strong preconfiguration on the 10 smallest Walshaw graphs

KAFFPA="$(dirname "$0")/build/kaffpa"
GRAPH_DIR="/home/c_schulz/projects/graph_collection/all_walschaw_graphs"
VERIFY="/tmp/verify_connected.py"
TMPDIR="/tmp/connected_blocks_test"
mkdir -p "$TMPDIR"

# Find 10 smallest graphs by file size
GRAPHS=$(ls -S "$GRAPH_DIR"/*.graph | tail -10 | tac)

K_VALUES="2 4 8 16 32 64"

PASS=0
FAIL=0
TOTAL=0

echo "============================================="
echo "Connected Blocks Test Suite (strong only)"
echo "============================================="
echo ""

for graph in $GRAPHS; do
    gname=$(basename "$graph" .graph)
    for k in $K_VALUES; do
        outfile="$TMPDIR/${gname}_strong_k${k}.txt"
        TOTAL=$((TOTAL + 1))

        # Run kaffpa
        result=$("$KAFFPA" "$graph" --k "$k" --preconfiguration strong \
            --connected_blocks --output_filename "$outfile" 2>&1)
        cut=$(echo "$result" | grep "^cut" | awk '{print $2}')
        balance=$(echo "$result" | grep "^balance" | awk '{print $2}')

        # Verify connectivity
        verify=$(python3 "$VERIFY" "$graph" "$outfile" 2>&1)
        if echo "$verify" | grep -q "^PASS"; then
            status="PASS"
            PASS=$((PASS + 1))
        else
            status="FAIL"
            FAIL=$((FAIL + 1))
        fi

        printf "%-12s k=%-3d cut=%-8s balance=%-8s %s\n" \
            "$gname" "$k" "$cut" "$balance" "$status"
    done
done

echo ""
echo "============================================="
echo "Results: $PASS/$TOTAL passed, $FAIL failed"
echo "============================================="

exit $FAIL
