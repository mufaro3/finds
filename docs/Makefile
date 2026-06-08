SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help clean latexpdf apidoc Makefile

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf $(BUILDDIR)
	rm -rf $(API_DIR)

%: Makefile apidoc
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

latexpdf:
	@echo "Building LaTeX..."
	@mkdir -p /tmp/sphinx-latex
	@python3 -m sphinx -b latex source /tmp/sphinx-latex
	@$(MAKE) -C /tmp/sphinx-latex all-pdf
	@cp /tmp/sphinx-latex/*.pdf build/
	@cp /tmp/sphinx-latex/finds.log build/
