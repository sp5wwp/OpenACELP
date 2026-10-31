# OpenACELP build
#
#   make                  build the codec library (build/libopenacelp.a) and unit tests
#   make apps             build oa_encode / oa_decode
#   make test             run the unit tests
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

.PHONY: all clean regen-tables apps test consts

all: $(LIB) $(BUILD)/test_unit

$(BUILD):
	mkdir -p $@

$(BUILD)/%.o: src/%.c src/oa_dsp.h src/oa_consts.h include/openacelp.h | $(BUILD)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(LIB): $(LIB_OBJ)
	$(AR) rcs $@ $^

# regenerate the formula-based constant tables (committed)
consts:
	python3 tools/gen_consts.py src/oa_consts.h

$(BUILD)/test_unit: tests/test_unit.c $(LIB)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIB) -lm

test: $(BUILD)/test_unit
	$(BUILD)/test_unit

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
	rm -rf $(BUILD)/*.o $(BUILD)/*.a $(BUILD)/oa_* $(BUILD)/test_unit
