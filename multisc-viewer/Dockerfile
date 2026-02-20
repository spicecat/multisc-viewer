FROM node:lts-bookworm-slim
LABEL org.opencontainers.image.title="multisc-viewer" \
	org.opencontainers.image.description="Web application for MultiSC-Viewer"
WORKDIR /app
COPY package.json pnpm-lock.yaml* ./
RUN npm install -g pnpm && pnpm install --frozen-lockfile
COPY . .
RUN pnpm run build && pnpm prune --prod
ENV NODE_ENV=production PORT=3000 MULTISC_VIEWER_DATA_DIR=/data
EXPOSE 3000
CMD ["node","build/index.js"]
