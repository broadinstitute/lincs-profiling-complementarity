import sys
import subprocess
import pkg_resources

required = {'iterative-stratification', 'pytorch-tabnet', 'umap-learn', 'scikit-learn', 'matplotlib', 'seaborn',
            'https://github.com/Phlya/adjustText/archive/master.zip','pandas','numpy', 'scikit-multilearn'} 
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed


if missing:
    # implement pip as a subprocess:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install',*missing])