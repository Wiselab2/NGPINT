# Checking recipes locally
## Fork the bioconda github 
git clone https://github.com/sagnikbanerjee15/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git

## Make sure our master is up to date with Bioconda
git checkout master
git pull upstream master
git push origin master

## Create and checkout a new branch for our work
git checkout -b update_my_recipe_from_ceres

cd recipes
conda skeleton pypi finder
Make changes in the finder/meta.yaml file and add build.sh
cd ..

./bootstrap.py --no-docker /project/maizegdb/sagnik/FINDER/miniconda
source ~/.config/bioconda/activate
bioconda-utils lint recipes config.yml --loglevel debug --full-report --git-range master HEAD
bioconda-utils build recipes config.yml --mulled-test --git-range master HEAD

# Choose the edited files to commit:
git add ngpint
# Create a commit (will open editor)
git commit -m "Push"
# Push your commit to GitHub
git push --set-upstream origin update_my_recipe_from_ceres

# Create a pull request
Now head back to GitHub. Go to your fork of the bioconda-recipes repository and find the branch you have just created in the Branch: drop down. You should now see a message saying This branch is 1 commit ahead [...] bioconda:master. To the right of that line you will find a button Pull Request. Click this and follow the instructions to open a new Pull Request.

# Push to my github

```bash
# Create a new repository on the command line
git init
git add --all
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/Wiselab2/NGPINT.git
git push -u origin main

# Push an existing repository from the command line
git remote add origin https://github.com/Wiselab2/NGPINT.git
git branch -M main
git push -u origin main
```



