# MultiSC-Viewer

Web application for [MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer). 

## Run locally

Requires [Node.js](https://nodejs.org/) and [PM2](https://pm2.keymetrics.io/).

```bash
# Clone the repository
git clone https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer/multisc-viewer

# Install dependencies
npm install

# Start the app
npm run dev
```

View app at <http://localhost:5173>.

## Docker

```bash
docker build -t multisc-viewer .
docker run -d -p 3000:80 multisc-viewer
```

## Datasets

## Publications