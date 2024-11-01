# CassetteCrafter

CassetteCrafter is a web app for generating cassettes for Golden Gate assembly. Users input DNA sequences, plasmid backbones, and protein sequences, and the app outputs optimized sequences ready for ordering.

## Project Structure

- `app_script.py`: Main Flask app script with routes and backend processing.
- `templates/`: 
  - `index.html`: Input page for user data.
  - `result.html`: Result page for optimized sequence download.
- `static/`: 
  - `ggpic.png`: Image displayed on the result page.

## Program Flow

1. **Input Page (`index.html`)**: Users input data for cassette generation.
2. **Backend Processing (`app_script.py`)**: Flask processes inputs and renders the result page with optimized sequence.
3. **Result Page (`result.html`)**: Displays the optimized sequence with a download link and the image `ggpic.png`.

## Installation and Setup

### Requirements

- Python 3.x
- Flask (`pip install flask`)

### Setup

1. Clone the repository:
   git clone https://github.com/yourusername/cassettecrafter.git
   cd cassettecrafter

2. Instal Flask
   pip install flask

3. Run the App
   python app_script.py

4. Open Web App
   Open a web browser and go to http://127.0.0.1:5000
