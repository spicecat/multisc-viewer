FROM node:current-alpine
LABEL org.opencontainers.image.title="multisc-viewer" \
	org.opencontainers.image.description="Web application for Multisc-Viewer"
WORKDIR /app
COPY package.json pnpm-lock.yaml* ./
RUN npm install -g pnpm && pnpm install --frozen-lockfile
COPY . .
RUN pnpm run build && pnpm prune --prod
ENV NODE_ENV=production PORT=3000 DATASETS_DIR=/data/datasets PUBLICATIONS_DIR=/data/publications PLOTS_DIR=/data/plots
EXPOSE 3000
CMD ["node","build/index.js"]
