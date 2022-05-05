env:
	# create and configurate conda environment
	conda env create -f environment.yml

html:
	jupyter book build .

html-hub:
	jupyter-book config sphinx .
	sphinx-build . _build/html -D html_baseurl=${JUPYTERHUB_SERVICE_PREFIX}/proxy/absolute/8000

clean:
	rm audio/*
	rm figures/*
	rm -rf _build