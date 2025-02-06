#!/bin/bash

# Remove annoying DS store files that pose cybersecurity risks.
#!/bin/bash

# Remove all .DS_Store files in the current directory and subdirectories
find . -type f -name ".DS_Store" -exec rm -f {} +

# Optional: Add .DS_Store to global .gitignore to prevent future tracking
if [ -f ~/.gitignore_global ]; then
    if ! grep -q ".DS_Store" ~/.gitignore_global; then
        echo ".DS_Store" >> ~/.gitignore_global
        git config --global core.excludesfile ~/.gitignore_global
        echo "Added .DS_Store to global .gitignore."
    fi
fi

echo "All .DS_Store files removed."
