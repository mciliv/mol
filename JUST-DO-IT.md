# JUST DO IT - Webkit Extraction Rule

## The Rule
When extracting webkit as a submodule, just do it. Don't overthink, don't create bureaucracy.

## Steps (in order, no exceptions)

1. **Create webkit repo** (if not exists)
   ```bash
   # Create new repo on GitHub/GitLab/etc
   # Name it: webkit
   # Make it public
   ```

2. **Clone webkit repo locally**
   ```bash
   git clone <webkit-repo-url> ../webkit
   ```

3. **Copy generic files to webkit**
   ```bash
   cp webkit.sh README-webkit.md ../webkit/
   cd ../webkit
   git add . && git commit -m "Add webkit files" && git push
   ```

4. **Remove files from mol project**
   ```bash
   cd ../mol
   rm webkit.sh README-webkit.md
   git add -A && git commit -m "Remove webkit files (now submodule)"
   ```

5. **Add webkit as submodule**
   ```bash
   git submodule add <webkit-repo-url> webkit
   git commit -m "Add webkit submodule"
   git push
   ```

6. **Update any scripts that use webkit**
   ```bash
   # Change: source webkit.sh
   # To: source webkit/webkit.sh
   ```

## That's it. No more steps.

## Usage
```bash
# In any script
source webkit/webkit.sh
log "Hello world"
gcloud_deploy "my-function" "us-central1" "./src"
```

## Remember
- Don't create complex structures
- Don't add unnecessary files
- Don't overthink the process
- Just do it 