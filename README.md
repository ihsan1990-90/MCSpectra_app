# How to Clone and Run MCSpectra locally
MCSpectra is online at www.mcspectra.com but can be cloned and run locally through the following steps:
## Step 1: Install Prerequisites
Before you start, ensure you have the following installed:
1.	**Python (3.8 or later)**
-	Download and install from [python.org](python.org)
-	Verify installation: `python --version`
2.	Git
-	Download and install from [git-scm.com](git-scm.com)
-	Verify installation: `git --version`
---
## Step 2: Clone the MCSpectra GitHub Repository
1.	Open **Terminal** (macOS/Linux) or **Command Prompt/PowerShell** (Windows).
2.	Navigate to the directory where you want to store the project: 
<pre> >> cd path/to/your/directory </pre>
3.	Clone the repository: 
<pre> >> git clone https://github.com/ihsan1990-90/MCSpectra_app </pre>
4.	Move into the project folder: 
<pre> >> cd MCSpectra_app </pre>
---
## Step 3: Set Up a Virtual Environment (Recommended)
To avoid dependency conflicts, create a virtual environment:
**Windows (Command Prompt/PowerShell)**
<pre> >> python -m venv venv </pre>
<pre> >> venv\Scripts\activate </pre>
**macOS/Linux (Terminal)**
<pre> >> python3 -m venv venv </pre>
<pre> >> source venv/bin/activate </pre>
---
## Step 4: Install Dependencies
Run the following command to install the required packages:
<pre> >> pip install -r requirements.txt </pre>
---
## Step 5: Run MCSpectra
Execute the following command inside the project directory:
<pre> >> streamlit run spectra_w_uncertainty_SL.py </pre>
---
## Step 6: Open the App in Browser
Once the app starts, you’ll see an output like:
<pre> You can now view your Streamlit app in your browser.
   Local URL: http://localhost:8501
   Network URL: http://192.168.x.x:8501 </pre>
Open the Local URL in a web browser to view the app. 
________________________________________

