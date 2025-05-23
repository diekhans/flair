# Releasing FLAIR bioconda package

Full details are here: https://bioconda.github.io/contributor/index.html

0. Create a pip package (this is used by conda)
1. Fork the bioconda recipes repo: https://github.com/bioconda/bioconda-recipes/fork
2. git clone that directory to your local computer
3. create a bioconda environment for testing: 
   - conda create -n bioconda -c conda-forge -c bioconda bioconda-utils pytorch
   - conda activate bioconda
4. update the recipe in recipes/flair/meta.yaml
5. validate the recipe using:
   - bioconda-utils lint --packages flair
6. test the recipe using:
   - bioconda-utils build --docker --mulled-test --force --packages flair >&bc-docker.log
6. test the recipe using
   - bioconda-utils build --mulled-test --packages flair 2>& bc-pkg.log
7. if no errors are found:
   -git commit/git push to your repo
8. submit a pull request via https://github.com/bioconda/bioconda-recipes/
   - This starts a testing process, please follow instructions on the page.
