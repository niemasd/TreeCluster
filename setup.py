from setuptools import setup,find_packages

setup(
        name='treecluster',    # This is the name of your PyPI-package.
        version='1.0.1',    # Update the version number for new releases
        scripts=['TreeCluster.py',], # The name of your script, and also the command you'll be using for calling it
        description='TreeCluster: a tool for clustering biological sequences using phylogenetic trees.',
        long_description='TreeCluster is a tool that, given a tree T (Newick format) and a distance threshold t, \
         finds the minimum number of clusters of the leaves of T such that some user-specified constraint is met \
         in each cluster. ',
        long_description_content_type='text/plain',
        url='https://github.com/niemasd/TreeCluster',
        author='Niema Moshiri',
        author_email='niemamoshiri@gmail.com',
        packages=find_packages(),
        zip_safe = False,
        install_requires=['niemads','treeswift'],
        include_package_data=True
)
