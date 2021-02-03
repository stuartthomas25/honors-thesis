pdf: notes.md assets/zotero.yaml
	#pandoc --filter pandoc-citeproc -o notes.pdf notes.md
	pandoc --resource-path=assets --default-image-extension=png --citeproc --from markdown+implicit_figures -o notes.pdf notes.md
