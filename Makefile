# OpenACELP build
#
#   make                  build the codec library (build/libopenacelp.a) and unit tests
#   make tetra-tables     fetch ETSI EN 300 395-2 and extract the TETRA tables
#   make tetra-tables ETSI_ZIP=path/to/en_30039502v010303p0.zip
#                         same, from a ZIP you already downloaded
#   make apps             build oa_encode / oa_decode (needs the TETRA tables)
#   make test             run the unit tests
#
# The TETRA tables are ETSI material. They are generated into $(BUILD)/ only,
# which is git-ignored. Never commit them, and never distribute a binary
# built with them.

CC      ?= gcc
CFLAGS  ?= -std=gnu17 -O2 -Wall -Wextra
BUILD   := build
INC     := -Iinclude -Isrc -I$(BUILD)

LIB_SRC := src/oa_dsp.c src/oa_bits.c src/oa_dec.c src/oa_enc.c
LIB_OBJ := $(LIB_SRC:src/%.c=$(BUILD)/%.o)
LIB     := $(BUILD)/libopenacelp.a

ETSI_VER    := 01.03.03
ETSI_TAG    := $(subst .,,$(ETSI_VER))
ETSI_FILE   := en_30039502v$(ETSI_TAG)p0.zip
ETSI_URL    := https://www.etsi.org/deliver/etsi_en/300300_300399/30039502/$(ETSI_VER)_60/$(ETSI_FILE)
ETSI_SHA256 := 06fa7611fc92d46a2c0fb173aa3e007dc80a0c3f03b94eb8ed181c9b7e5f3034
ETSI_DIR    := $(BUILD)/etsi
ETSI_ZIP    ?= $(ETSI_DIR)/$(ETSI_FILE)

.PHONY: all clean tetra-tables apps test consts

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

# --- applications (need the TETRA tables) ---------------------------------

apps: $(BUILD)/oa_decode $(BUILD)/oa_encode

$(BUILD)/tetra_tables.o: $(BUILD)/tetra_tables.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(BUILD)/oa_%: apps/oa_%.c $(LIB) $(BUILD)/tetra_tables.o
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(BUILD)/tetra_tables.o $(LIB)

# --- TETRA tables (opt-in) -------------------------------------------------

tetra-tables: $(BUILD)/tetra_tables.c

$(ETSI_DIR)/$(ETSI_FILE):
	@echo "Downloading ETSI EN 300 395-2 $(ETSI_VER) from etsi.org."
	@echo "By continuing you accept ETSI's terms of use for this material."
	mkdir -p $(ETSI_DIR)
	@rm -f $@.part
	@if command -v wget >/dev/null 2>&1; then \
	  wget -O $@.part $(ETSI_URL); \
	else \
	  curl -fL --retry 2 -o $@.part $(ETSI_URL); \
	fi || { \
	  rm -f $@.part; \
	  echo ""; \
	  echo "Automatic download failed."; \
	  echo "Download the file in a web browser from:"; \
	  echo "  $(ETSI_URL)"; \
	  echo "then run:"; \
	  echo "  make tetra-tables ETSI_ZIP=/path/to/$(ETSI_FILE)"; \
	  exit 1; }
	mv $@.part $@

$(ETSI_DIR)/.verified: $(ETSI_ZIP)
	@mkdir -p $(ETSI_DIR)
	echo "$(ETSI_SHA256)  $<" | sha256sum -c -
	@touch $@

$(BUILD)/tetra_tables.c: $(ETSI_DIR)/.verified tools/extract_etsi_tables.py
	rm -rf $(ETSI_DIR)/src
	unzip -L -C -o -q $(ETSI_ZIP) 'c-code/*' -d $(ETSI_DIR)/src
	python3 tools/extract_etsi_tables.py $(ETSI_DIR)/src/c-code $(BUILD)

clean:
	rm -rf $(BUILD)/*.o $(BUILD)/*.a $(BUILD)/oa_* $(BUILD)/test_unit
