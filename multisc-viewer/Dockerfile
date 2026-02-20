FROM node:lts-bookworm-slim AS builder

WORKDIR /app

RUN corepack enable && corepack prepare pnpm@latest --activate
COPY package.json pnpm-lock.yaml ./
RUN pnpm install --frozen-lockfile --ignore-scripts

COPY . .

RUN pnpm rebuild
RUN pnpm run build
RUN pnpm prune --prod --ignore-scripts


FROM node:lts-bookworm-slim AS runner

LABEL org.opencontainers.image.title="multisc-viewer" \
	org.opencontainers.image.description="Web application for MultiSC-Viewer"

WORKDIR /app

ENV NODE_ENV=production \
	PORT=3000 \
	MULTISC_VIEWER_DATA_DIR=/data

COPY --from=builder /app/build ./build
COPY --from=builder /app/node_modules ./node_modules
COPY --from=builder /app/package.json ./package.json

EXPOSE 3000

CMD ["node", "build/index.js"]