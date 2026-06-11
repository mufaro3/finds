SPHINXOPTS  ?=
SPHINXBUILD ?= sphinx-build
SOURCEDIR   = source
BUILDDIR    = build

TMP_HTML    = /tmp/sphinx-html
TMP_LATEX   = /tmp/sphinx-latex

.PHONY: help clean html latexpdf all

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" \
	$(SPHINXOPTS) $(O)

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(TMP_HTML)
	rm -rf $(TMP_LATEX)

# ---------- HTML ----------
html:
	@echo "Building HTML docs..."
	@rm -rf $(TMP_HTML)
	@mkdir -p $(TMP_HTML)

	@python3 -m sphinx \
		-b html \
		"$(SOURCEDIR)" \
		"$(TMP_HTML)"

	@mkdir -p "$(BUILDDIR)/html"
	@cp -r $(TMP_HTML)/* "$(BUILDDIR)/html/"

# ---------- PDF ----------
latexpdf:
	@echo "Building PDF docs..."
	@rm -rf $(TMP_LATEX)
	@mkdir -p $(TMP_LATEX)

	@python3 -m sphinx \
		-b latex \
		"$(SOURCEDIR)" \
		"$(TMP_LATEX)"

	@$(MAKE) -C "$(TMP_LATEX)" all-pdf

	@mkdir -p "$(BUILDDIR)"
	@cp $(TMP_LATEX)/*.pdf "$(BUILDDIR)/"
	@cp $(TMP_LATEX)/*.log "$(BUILDDIR)/" || true

# ---------- Build everything ----------
all: html latexpdf
	@echo ""
	@echo "Build complete:"
	@echo "  HTML: $(BUILDDIR)/html/index.html"
	@echo "  PDF : $(BUILDDIR)/*.pdf"
