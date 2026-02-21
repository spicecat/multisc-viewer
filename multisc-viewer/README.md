# MultiSC-Viewer Web App

Web application for [MultiSC-Viewer](https://github.com/spicecat/multisc-viewer).

## Run Locally

Requires [Node.js](https://nodejs.org/).

```bash
# Clone the repository
git clone --single-branch https://github.com/spicecat/multisc-viewer.git
cd multisc-viewer/multisc-viewer

# Install dependencies
npm install

# Start the app
npm run dev
```

View app at <http://localhost:5173>.

## Docker

Requires [Docker](https://www.docker.com/).

```bash
docker build -t multisc-viewer .
docker run -d -p 3000:80 multisc-viewer
```

## Environment Variables

- `MULTISC_VIEWER_DATA_DIR` - Path to the data directory (default: `/data`)

---

## Related

- MultiSC: [https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer](https://github.com/spicecat/multisc-viewer/tree/main/README.md)
- Add MultiSC-Data: [https://github.com/spicecat/multisc-viewer/tree/main/data](https://github.com/spicecat/multisc-viewer/tree/main/data/README.md)
- Run MultiSC-Viewer Daemon: [https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR](https://github.com/spicecat/multisc-viewer/tree/main/MultiSCViewerR/README.md)
- Run MultiSC-Viewer: [https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer](https://github.com/spicecat/multisc-viewer/tree/main/multisc-viewer/README.md)
