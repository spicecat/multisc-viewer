# MultiSC-Viewer

Web application for [MultiSC-Viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main).

## Run Locally

Requires [Node.js](https://nodejs.org/).

```bash
# Clone the repository
git clone --single-branch -b main https://git.jasonxu.dev/JasonXu/plot-viewer.git
cd plot-viewer/multisc-viewer

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

---

## Related

- MultiSC: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/README.md)
- Add MultiSC-Data: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/data](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/data/README.md)
- Run MultiSC-Daemon: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-daemon](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-daemon/README.md)
- Run MultiSC-Viewer: [https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer](https://git.jasonxu.dev/JasonXu/plot-viewer/src/branch/main/multisc-viewer/README.md)
