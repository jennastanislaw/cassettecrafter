# from flask import Flask, render_template, request, send_file, session, redirect, url_for
# import io
# import os

# app = Flask(__name__)  # Initialize the Flask application
# app.secret_key = os.urandom(24)  # Secret key for session management, used to securely sign the session cookie

# def manipulate_sequences(backbone, protein):
#     # Ensure the sequences are exactly 10 base pairs long
#     if len(backbone) != 10 or len(protein) != 10:
#         raise ValueError("Both Backbone and Protein sequences must be exactly 10 base pairs long.")
    
#     # Create a summary of both sequences
#     summary = f"Backbone Sequence: {backbone}\nProtein Sequence: {protein}"
#     return summary

# @app.route('/')
# def home():
#     return render_template('index.html')  # Render the home page with the form

# @app.route('/optimize', methods=['POST'])
# def optimize():
#     backbone = request.form['backbone']  # Get the backbone sequence from the form
#     protein = request.form['protein']  # Get the protein sequence from the form
    
#     if not backbone or not protein:
#         return "Error: Both Backbone and Protein sequences are required", 400  # Return an error if either sequence is missing
    
#     try:
#         # Call the function to manipulate the sequences
#         result = manipulate_sequences(backbone, protein)
        
#         # Store the result in the session
#         session['result'] = result
#         return redirect(url_for('result'))  # Redirect to the result page
#     except ValueError as ve:
#         return str(ve), 400  # Return an error if the sequences are not 10 base pairs long
#     except Exception as e:
#         return f"An error occurred: {str(e)}", 500  # Return a generic error for any other exceptions

# @app.route('/result')
# def result():
#     if 'result' not in session:
#         return redirect(url_for('home'))  # Redirect to the home page if there is no result in the session
#     return render_template('result.html')  # Render the result page with the download button

# @app.route('/download')
# def download():
#     if 'result' not in session:
#         return redirect(url_for('home'))  # Redirect to the home page if there is no result in the session
    
#     result = session['result']  # Get the result from the session
#     output = io.StringIO(result)  # Create a file-like object with the result
#     return send_file(io.BytesIO(output.getvalue().encode()), 
#                      download_name="optimized_sequences.txt", 
#                      as_attachment=True)  # Send the file as a download

# if __name__ == '__main__':
#     app.run(debug=True)  # Run the Flask application in debug mode



# from flask import Flask, render_template, request, send_file, session, redirect, url_for
# import io
# import os

# app = Flask(__name__)  # Initialize the Flask application
# app.secret_key = os.urandom(24)  # Secret key for session management, used to securely sign the session cookie

# def manipulate_sequences(backbone, protein):
#     # Ensure the sequences are exactly 10 base pairs long
#     if len(backbone) != 10 or len(protein) != 10:
#         raise ValueError("Both Backbone and Protein sequences must be exactly 10 base pairs long.")
    
#     # Create a summary of both sequences
#     summary = f"Backbone Sequence: {backbone}\nProtein Sequence: {protein}"
#     return summary

# @app.route('/')
# def home():
#     return render_template('index.html')  # Render the home page with the form

# @app.route('/optimize', methods=['POST'])
# def optimize():
#     backbone = request.form['backbone']  # Get the backbone sequence from the form
#     protein = request.form['protein']  # Get the protein sequence from the form
    
#     if not backbone or not protein:
#         return "Error: Both Backbone and Protein sequences are required", 400  # Return an error if either sequence is missing
    
#     try:
#         # Call the function to manipulate the sequences
#         result = manipulate_sequences(backbone, protein)
        
#         # Store the result in the session
#         session['result'] = result
#         return redirect(url_for('result'))  # Redirect to the result page
#     except ValueError as ve:
#         return str(ve), 400  # Return an error if the sequences are not 10 base pairs long
#     except Exception as e:
#         return f"An error occurred: {str(e)}", 500  # Return a generic error for any other exceptions

# @app.route('/result')
# def result():
#     result = session.get('result', None)
#     if result is None:
#         return redirect(url_for('home'))
#     return f"<pre>{result}</pre>"

# if __name__ == '__main__':
#     app.run(debug=True)

from flask import Flask, render_template, request, send_file, session, redirect, url_for
import io
import os

app = Flask(__name__)  # Initialize the Flask application
app.secret_key = os.urandom(24)  # Secret key for session management, used to securely sign the session cookie

def manipulate_sequences(backbone, protein):
    # Ensure the sequences are exactly 10 base pairs long
    if len(backbone) != 10 or len(protein) != 10:
        raise ValueError("Both Backbone and Protein sequences must be exactly 10 base pairs long.")
    
    # Create a summary of both sequences
    summary = f"Backbone Sequence: {backbone}\nProtein Sequence: {protein}"
    return summary

@app.route('/')
def home():
    return render_template('index.html')  # Render the home page with the form

@app.route('/optimize', methods=['POST'])
def optimize():
    backbone = request.form['backbone']  # Get the backbone sequence from the form
    protein = request.form['protein']  # Get the protein sequence from the form
    
    if not backbone or not protein:
        return "Error: Both Backbone and Protein sequences are required", 400  # Return an error if either sequence is missing
    
    try:
        # Call the function to manipulate the sequences
        result = manipulate_sequences(backbone, protein)
        
        # Store the result in the session
        session['result'] = result
        return redirect(url_for('result'))  # Redirect to the result page
    except ValueError as ve:
        return str(ve), 400  # Return an error if the sequences are not 10 base pairs long
    except Exception as e:
        return f"An error occurred: {str(e)}", 500  # Return a generic error for any other exceptions

@app.route('/result')
def result():
    result = session.get('result', None)
    if result is None:
        return redirect(url_for('home'))
    return render_template('result.html', result=result)

@app.route('/download')
def download():
    result = session.get('result', None)
    if result is None:
        return redirect(url_for('home'))
    
    # Create a file-like object to send as a downloadable file
    output = io.StringIO(result)
    output.seek(0)
    
    return send_file(output, mimetype='text/plain', as_attachment=True, attachment_filename='result.txt')

if __name__ == '__main__':
    app.run(debug=True)
