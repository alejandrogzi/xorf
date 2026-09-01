#!/usr/bin/env bash
# ==============================================================================
# test_containers.sh - Test all ORF pipeline containers
# ==============================================================================
set -euo pipefail

# Configuration
REGISTRY="${DOCKER_REGISTRY:-localhost}"
TAG="${DOCKER_TAG:-latest}"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

PASSED=0
FAILED=0

# Test result tracking
test_result() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}"
        ((PASSED++)) || true
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}"
        ((FAILED++)) || true
        return 1
    fi
}

echo "=============================================="
echo "Testing ORF Pipeline Containers"
echo "=============================================="
echo ""

# ==============================================================================
# Test 1: orf-chunk
# ==============================================================================
echo -e "${BLUE}[TEST 1/12]${NC} orf-chunk: Help command"
echo "Running: docker run --rm "${REGISTRY}/orf-chunk:${TAG}" chunk --help"
docker run --rm "${REGISTRY}/orf-chunk:${TAG}" chunk --help > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 2/12]${NC} orf-chunk: Version check"
echo "Running: docker run --rm "${REGISTRY}/orf-chunk:${TAG}" --version"
docker run --rm "${REGISTRY}/orf-chunk:${TAG}" --version > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 3/12]${NC} orf-chunk: Binary size check"
echo "Running: docker run --rm "${REGISTRY}/orf-chunk:${TAG}" sh -c "du -h /usr/local/bin/orf | cut -f1""
SIZE=$(docker run --rm --entrypoint sh "${REGISTRY}/orf-chunk:${TAG}" -c "du -h /usr/local/bin/orf | cut -f1")
echo "  Binary size: ${SIZE}"
test_result

# ==============================================================================
# Test 2: orf-tai
# ==============================================================================
echo -e "${BLUE}[TEST 4/12]${NC} orf-tai: Help command"
echo "Running: docker run --rm "${REGISTRY}/orf-tai:${TAG}" tai --help"
docker run --rm "${REGISTRY}/orf-tai:${TAG}" tai --help > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 5/12]${NC} orf-tai: translationai availability"
echo "Running: docker run --rm --entrypoint sh "${REGISTRY}/orf-tai:${TAG}" -c "translationai --help""
docker run --rm --entrypoint sh "${REGISTRY}/orf-tai:${TAG}" -c "translationai --help" > /dev/null 2>&1
test_result

# echo -e "${BLUE}[TEST 6/12]${NC} orf-tai: Python environment check"
# echo "Running: docker run --rm --entrypoint sh "${REGISTRY}/orf-tai:${TAG}" -c "micromamba run -n taienv python --version""
# docker run --rm --entrypoint sh "${REGISTRY}/orf-tai:${TAG}" -c "micromamba run -n taienv python --version" 2>&1 | grep "Python 3.9"
# test_result

# ==============================================================================
# Test 3: orf-samba
# ==============================================================================
echo -e "${BLUE}[TEST 7/12]${NC} orf-samba: Help command"
echo "Running: docker run --rm "${REGISTRY}/orf-samba:${TAG}" samba --help"
docker run --rm "${REGISTRY}/orf-samba:${TAG}" samba --help > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 8/12]${NC} orf-samba: rnasamba availability"
echo "Running: docker run --rm --entrypoint "${REGISTRY}/orf-samba:${TAG}" -c "rnasamba --help""
docker run --rm --entrypoint sh "${REGISTRY}/orf-samba:${TAG}" -c "rnasamba --help" > /dev/null 2>&1
test_result

# echo -e "${BLUE}[TEST 9/12]${NC} orf-samba: Python environment check"
# echo "Running: docker run --rm  --entrypoint sh "${REGISTRY}/orf-samba:${TAG}" -c "micromamba run -n sambaenv python --version""
# docker run --rm  --entrypoint sh "${REGISTRY}/orf-samba:${TAG}" -c "micromamba run -n sambaenv python --version" 2>&1 | grep "Python 3.6"
# test_result

# ==============================================================================
# Test 4: orf-blast
# ==============================================================================
echo -e "${BLUE}[TEST 10/12]${NC} orf-blast: Help command"
echo "Running: docker run --rm "${REGISTRY}/orf-blast:${TAG}" blast --help"
docker run --rm "${REGISTRY}/orf-blast:${TAG}" blast --help > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 11/12]${NC} orf-blast: diamond availability"
echo "Running: docker run --rm  --entrypoint sh "${REGISTRY}/orf-blast:${TAG}" -c "diamond --version""
docker run --rm  --entrypoint sh "${REGISTRY}/orf-blast:${TAG}" -c "diamond --version" > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 12/13]${NC} orf-blast: orfipy availability"
echo "Running: docker run --rm --entrypoint sh "${REGISTRY}/orf-blast:${TAG}" -c "orfipy --version""
docker run --rm --entrypoint sh "${REGISTRY}/orf-blast:${TAG}" -c "orfipy --version" > /dev/null 2>&1
test_result

echo -e "${BLUE}[TEST 13/13]${NC} orf-blast: mmseqs availability"
echo "Running: docker run --rm --entrypoint sh "${REGISTRY}/orf-blast:${TAG}" -c "mmseqs version""
docker run --rm --entrypoint sh "${REGISTRY}/orf-blast:${TAG}" -c "mmseqs version" > /dev/null 2>&1
test_result

# ==============================================================================
# Summary
# ==============================================================================
echo ""
echo "=============================================="
echo "Test Summary"
echo "=============================================="
echo -e "Total tests: $((PASSED + FAILED))"
echo -e "${GREEN}Passed: ${PASSED}${NC}"
echo -e "${RED}Failed: ${FAILED}${NC}"
echo ""

if [ ${FAILED} -eq 0 ]; then
    echo -e "${GREEN}All tests passed! ✓${NC}"
    echo ""
    echo "Container sizes:"
    docker images | grep "orf-" | grep "${TAG}" | awk '{print $1 ":" $2 " - " $7 $8}'
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
