# Molecular Reality (mol)

This project includes a Node.js server that hosts a static front-end (HTML, CSS, JavaScript) and provides back-end API endpoints for processing molecular data using the OpenAI API, as well as Python utilities for SDF file generation.

## Prerequisites

- [Node.js](https://nodejs.org/) (>= 18)
- [Python](https://www.python.org/) (>= 3.8)
- [poetry](https://python-poetry.org/) (for Python dependencies)

## Python Environment Setup

To set up the Python virtual environment and install dependencies:

```bash
chmod +x scripts/helper.sh
./scripts/helper.sh
```

## Node.js Setup

Install Node.js dependencies:

```bash
npm install
```

### Environment Variables

Set your OpenAI API key:

```bash
export OPENAI_API_KEY=<your_api_key_here>
```

### Running the Server

Start the Node.js server (serves both front-end and back-end API):

```bash
npm start
```

For development with automatic reloads:

```bash
npm run dev
```

The application will be available at [http://localhost:3000](http://localhost:3000).

## API Endpoints

- `POST /list-molecules`: Analyze an image click and return possible molecule structures.
- `POST /list-molecules-text`: Provide an object description and get molecule SMILES in JSON.
