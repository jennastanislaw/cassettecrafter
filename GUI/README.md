# CassetteCrafter

CassetteCrafter is a web app for generating plasmid cassettes for Golden Gate assembly. Users input DNA sequences, plasmid backbones, and protein sequences, and the app outputs optimized sequences ready for ordering.

## Project Structure

- `app_script.py`: Main Flask app script with routes and backend processing.
- `temp/`:
   - Data files
- `templates/`: 
  - `index.html`: Input page for user data.
  - `gg_assembly.html`: Information page.
  - `result.html`: Result page for optimized sequence download.
- `static/`: 
  - Images used in the html files.
- `Tests`: 
   - `test_app_script.py`: Result page for optimized sequence download.

## Program Flow

1. **Input Page (`index.html`)**: Users input data for cassette generation.
2. **Information Page (`gg_assembly.html`)**: Users can learn about the golden gate assembly process.
3. **Backend Processing (`app_script.py`)**: Flask processes inputs and renders the result page with optimized sequence.
4. **Result Page (`result.html`)**: Displays the optimized sequence with a download link 

## Installation and Setup

### Requirements

- Python 3.x
- Flask (`pip install flask`)
- biopython
- pandas

### Setup

1. Clone the repository:
   git clone https://github.com/yourusername/cassettecrafter.git
   cd cassettecrafter

2. Instal Flask
   pip install flask

3. Run the App 
   Ensure the app is run from the cassettecrafter directory
   python3 GUI/app_script.py

4. Open Web App
   Open a web browser and go to http://127.0.0.1:5000
