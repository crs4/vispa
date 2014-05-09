EXPORT_DIR = svn_export
GENERATED_FILES = AUTHORS MANIFEST README bl/tiget/version.py
PY_V := $(shell python -c 'import sys; print "%d.%d" % sys.version_info[:2]')

.PHONY: all build build_py install install_py install_user install_user_py dist clean distclean uninstall_user

all: build

build:
	python setup.py build

build_py:
	python setup.py build_py

install: build
	python setup.py install --skip-build

install_py: build_py
	python setup.py install --skip-build

install_user: build
	python setup.py install --skip-build --user

install_user_py: build_py
	python setup.py install --skip-build --user

dist:
	rm -rf $(EXPORT_DIR) && svn export . $(EXPORT_DIR)
	cd $(EXPORT_DIR) && python setup.py sdist -k

clean:
	rm -rf build
	rm -f $(GENERATED_FILES)
	find . -regex '.*\(\.pyc\|\.pyo\|~\|\.so\)' -exec rm -fv {} \;

distclean: clean
	rm -rf $(EXPORT_DIR) dist

uninstall_user:
	rm -rf ~/.local/lib/python$(PY_V)/site-packages/bl/tiget
	rm -f ~/.local/lib/python$(PY_V)/site-packages/tiget*
