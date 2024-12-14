'''
This module generates  Graphical User Interface (GUI) for the assembly library generation. 
It uses Flask to create a web application that allows users to upload gene and mutation files,
and then generates a Golden Gate assembly library based on the uploaded files.

The main functions of this module are:
    - home: Renders the home page with the form for uploading files and setting parameters
    - optimize: Handles the form submission and calls the generate_assembly_library function
    - result: Renders the result page with the generated library
    - download: Allows users to download the generated library result as a CSV file

Additional routes are defined for the Golden Gate assembly information page (gg_assmebly) 
and for the main function.
'''

from flask import Flask, render_template, request, send_file, session, redirect, url_for
import os
import io
import sys

# Add the src directory to the Python path
src_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../src"))

print(f"Adding src directory to path: {src_dir}")
sys.path.append(src_dir)
try:
    from main import generate_assembly_library  # Import the function from main.py
except ImportError as e:
    raise(e)
print("All modules imported successfully!")

app = Flask(__name__)
app.secret_key = os.urandom(24)

# Define the routes
@app.route('/')
def home():
    return render_template('index.html')  # Render the home page with the form

# Define the optimize route
@app.route('/optimize', methods=['POST'])
def optimize():
    # Get form inputs and files
    gene_file = request.files.get('gene_file')
    mutations_file = request.files.get('mutations_file')
    enzyme_name = request.form.get('enzyme_name', 'BbsI')  
    min_oligo_size = request.form.get('min_oligo_size')
    min_oligo_size = int(min_oligo_size) if min_oligo_size else 20
    max_oligo_size = int(request.form.get('max_oligo_size', 100))

    # Check if both files are uploaded
    if not gene_file or not mutations_file:
        return "Error: Both gene file and mutations file are required", 400
    
    try:
        # Save files temporarily
        temp_dir = './temp'
        os.makedirs(temp_dir, exist_ok=True)
        gene_file_path = os.path.join(temp_dir, gene_file.filename)
        mutations_file_path = os.path.join(temp_dir, mutations_file.filename)
        gene_file.save(gene_file_path)
        mutations_file.save(mutations_file_path)

        # Call the function generate_assembly_library from main.py
        final_df = generate_assembly_library(
            gene_file_path, mutations_file_path, enzyme_name, min_oligo_size, max_oligo_size
        )

        # Convert the result to CSV for download
        output = io.StringIO()
        final_df.to_csv(output, index=False)
        output.seek(0)
        # Save result to session
        session['result_csv'] = output.getvalue()  
        # Redirect to the result page
        return redirect(url_for('result'))
    except Exception as e:
        return f"An error occurred: {str(e)}", 500
    
        
# Define the gg_assembly route for the Golden Gate assembly info page
@app.route('/gg_assembly')
def gg_assembly():
    return render_template('gg_assembly.html')

# Define the result route to display the result page
@app.route('/result')
def result():
    result_csv = session.get('result_csv', None)
    if not result_csv:
        return redirect(url_for('home'))
    return render_template('result.html')  # Display the result page

# Define the download route to download the generated librar result as a CSV file
@app.route('/download')
def download():
    result_csv = session.get('result_csv', None)
    if not result_csv:
        return redirect(url_for('home'))

    # Convert string to bytes and wrap it with BytesIO
    output = io.BytesIO(result_csv.encode('utf-8'))  # Encode the string as bytes
    output.seek(0)  # Rewind the stream to the beginning

    # Send the file as a download
    return send_file(output, mimetype='text/csv', as_attachment=True, download_name='result.csv')

# Run the app
if __name__ == '__main__':
    app.run(debug=True)


