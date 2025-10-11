# Research Lab Book
# Lab Book Entry — 2025-10-07 
**git-tag:** 2.47.0.windows.1
 

---

### **Objective**
Set up the initial research project structure, create a digital lab book, initialize version control, and connect the repository to the institutional GitLab server (`gitlab.rlp.net`).

---

### **Work Steps**

**1. Created local project folder and structure**
```bash
mkdir -p Docs/{Material,Methods,Results}
```
Added initial lab book:
```bash
touch Docs/LabBook.md
```

---

**2. Initialized Git repository and first commit**
```bash
git init
git add .
git commit -m "Initial commit with folder structure and lab book"
```

---

**3. Configured Git remote connection to institutional GitLab**
Set the correct remote URL:
```bash
git remote set-url origin https://gitlab.rlp.net/a.hallab/tensor-omics.git 
```
Checked the connection:
```bash
git remote -v
```
Output:
```
origin  https://gitlab.rlp.net/a.hallab/tensor-omics.git (fetch)
origin  https://gitlab.rlp.net/a.hallab/tensor-omics.git (push)

```

---

**4. Pushed initial project to GitLab**
```bash
git push 
```
Verified upload by opening repository in browser:  
🔗 https://gitlab.rlp.net/a.hallab/tensor-omics/-/tree/feature/initial_submodules

---



### **Results**
- Project successfully uploaded to GitLab.  
- `Docs/LabBook.md` contains the first lab entry.  
- All project folders are visible on the GitLab web interface.  
- Git remote correctly configured for institutional access.

---

### **Next Steps**
- Start documenting daily research tasks inside `Docs/LabBook.md`.  
- Prepare first experiment scripts inside `Methods/`.  

---

### **Problems / Notes**
- Initially received `Could not resolve host: git.th-bingen.de` error due to wrong remote URL.  
- Resolved by replacing with the correct GitLab domain (`gitlab.rlp.net`).  
- Verified repository access and successful push.

---

### **Time Spent**
~3.0 hours (project setup, troubleshooting, documentation)
