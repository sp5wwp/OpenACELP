# OpenACELP build
#
#   make                  build the codec library (build/libopenacelp.a) and unit tests
#   make apps             build oa_encode / oa_decode
#   make test             run the unit tests and the decoder regression test
#   make conformance ETSI_TEST_DIR=dir
#                         EN 300 395-4 decoder test: TE1_46.CHD must decode to
#                         TE2_128.PCM byte for byte
#   make cortex-m4        compile the library for Cortex-M4 (CROSS=arm-none-eabi-)
#   make regen-tables ETSI_ZIP=path/to/en_30039502v010303p0.zip
#                         maintainers: regenerate src/oa_tables_tetra.c from
#                         the ETSI ZIP (to check it or to change its format)

CC      ?= gcc
CFLAGS  ?= -std=gnu17 -O2 -Wall -Wextra
BUILD   := build
INC     := -Iinclude -Isrc -I$(BUILD)

LIB_SRC := src/oa_dsp.c src/oa_bits.c src/oa_dec.c src/oa_enc.c src/oa_tables_tetra.c
LIB_OBJ := $(LIB_SRC:src/%.c=$(BUILD)/%.o)
LIB     := $(BUILD)/libopenacelp.a

ETSI_SHA256 := 06fa7611fc92d46a2c0fb173aa3e007dc80a0c3f03b94eb8ed181c9b7e5f3034
ETSI_DIR    := $(BUILD)/etsi

.PHONY: all clean regen-tables apps test consts conformance cortex-m4

all: $(LIB) $(BUILD)/test_unit $(BUILD)/test_tetra

$(BUILD):
	mkdir -p $@

HDRS := src/oa_dsp.h src/oa_consts.h src/oa_basicop.h src/oa_xdec.h include/openacelp.h

$(BUILD)/%.o: src/%.c $(HDRS) | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(LIB): $(LIB_OBJ)
	$(AR) rcs $@ $^

# regenerate the formula-based constant tables (committed)
consts:
	python3 tools/gen_consts.py src/oa_consts.h

$(BUILD)/test_unit: tests/test_unit.c $(LIB)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIB) -lm

$(BUILD)/test_tetra: tests/test_tetra.c tests/test_tetra_hashes.inc $(LIB)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIB)

test: $(BUILD)/test_unit $(BUILD)/test_tetra
	$(BUILD)/test_unit
	$(BUILD)/test_tetra

# --- EN 300 395-4 decoder conformance (test files not distributed) ---------

conformance: $(BUILD)/oa_decode
	@test -n "$(ETSI_TEST_DIR)" || { echo "usage: make conformance ETSI_TEST_DIR=dir-with-TE1_46.CHD-and-TE2_128.PCM"; exit 1; }
	$(BUILD)/oa_decode $(ETSI_TEST_DIR)/TE1_46.CHD $(BUILD)/TE1_46.PCM
	cmp $(BUILD)/TE1_46.PCM $(ETSI_TEST_DIR)/TE2_128.PCM && echo "conformance: TE1_46.CHD -> TE2_128.PCM bit-exact"

# --- Cortex-M4 compile check -------------------------------------------------

CROSS   ?= arm-none-eabi-
M4FLAGS ?= -mcpu=cortex-m4 -mthumb -mfloat-abi=soft -O2 -std=gnu17 -Wall -Wextra -ffunction-sections -fdata-sections

cortex-m4: | $(BUILD)
	mkdir -p $(BUILD)/m4
	for f in $(LIB_SRC); do $(CROSS)gcc $(M4FLAGS) $(INC) -c -o $(BUILD)/m4/$$(basename $$f .c).o $$f || exit 1; done
	$(CROSS)size $(BUILD)/m4/*.o

# --- applications -----------------------------------------------------------

apps: $(BUILD)/oa_decode $(BUILD)/oa_encode

$(BUILD)/oa_%: apps/oa_%.c $(LIB)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIB)

# --- maintainers: regenerate src/oa_tables_tetra.c ---------------------------

regen-tables:
	@test -n "$(ETSI_ZIP)" || { echo "usage: make regen-tables ETSI_ZIP=path/to/en_30039502v010303p0.zip"; exit 1; }
	echo "$(ETSI_SHA256)  $(ETSI_ZIP)" | sha256sum -c -
	rm -rf $(ETSI_DIR) && mkdir -p $(ETSI_DIR)
	unzip -L -C -o -q $(ETSI_ZIP) 'c-code/*' -d $(ETSI_DIR)
	python3 tools/extract_etsi_tables.py $(ETSI_DIR)/c-code src/oa_tables_tetra.c

clean:
	rm -rf $(BUILD)/*.o $(BUILD)/*.a $(BUILD)/oa_* $(BUILD)/test_unit $(BUILD)/test_tetra $(BUILD)/m4
