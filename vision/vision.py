from flask import Flask, request, jsonify
from flask_cors import CORS
from dotenv import load_dotenv
import os
import requests
import base64

# Load environment variables
load_dotenv()

app = Flask(__name__)
CORS(app)

# API Keys
OPENAI_API_KEY = 'test'
IMGBB_API_KEY = 'test' # Replace with your ImgBB API key

# os.getenv('OPENAI_API_KEY')
# if not OPENAI_API_KEY:
#     raise ValueError("OpenAI API key not found in environment variables")

@app.route('/analyze-image', methods=['POST'])
def analyze_image():
    try:
        # Get image data from request
        data = request.json
        if not data or 'image' not in data:
            return jsonify({'error': 'No image data provided'}), 400


        # Upload image to ImgBB
        imgbb_url = 'https://api.imgbb.com/1/upload'
        imgbb_payload = {
            'key': IMGBB_API_KEY,
            'image': data['image']
        }

        imgbb_response = requests.post(imgbb_url, data=imgbb_payload)
        if not imgbb_response.ok:
            return jsonify({'error': 'Failed to upload image'}), 500

        image_url = imgbb_response.json()['data']['url']

        # Prepare the request to OpenAI API
        headers = {
            'Authorization': f'Bearer {OPENAI_API_KEY}',
            'Content-Type': 'application/json'
        }

        payload = {
            'model': 'gpt-4',
            'messages': [
                {
                    'role': 'user',
                    'content': [
                        {
                            'type': 'text',
                            'text': f'What do you see in this image? You are a chemist helping us for education purposes. I am curious about the materials around us - from everyday objects to complex structures. I want to know what everything is made out of. Help me. First identify one main object in the image (around these coordinates: X: {data.get("coordinates", {}).get("x")}, Y: {data.get("coordinates", {}).get("y")}.) and then analyze any chemical compounds and answer with the compounds listed as SMILES in a json array dont have json markdown-- nothing else, dont provide any other info.X:0 Y:0 is the top left corner of the image. The image is 1000x1000 pixels.'
                        },
                        {
                            'type': 'image_url',
                            'image_url': {
                                'url': image_url
                            }
                        }
                    ]
                }
            ],
            'max_tokens': 300
        }

        # Make request to OpenAI API
        response = requests.post(
            'https://api.openai.com/v1/chat/completions',
            headers=headers,
            json=payload
        )

        # Check if the request was successful
        response.raise_for_status()
        result = response.json()

        print(result['choices'][0]['message']['content'])

        return jsonify({
            'analysis': result['choices'][0]['message']['content']
        })

    except requests.exceptions.RequestException as e:
        print(e)
        return jsonify({'error': f'Error calling OpenAI API: {str(e)}'}), 500
    except Exception as e:
        print(e)
        return jsonify({'error': f'Server error: {str(e)}'}), 500

if __name__ == '__main__':
    # For development, use adhoc certificates
    app.run(host='0.0.0.0', port=5001, debug=True)
