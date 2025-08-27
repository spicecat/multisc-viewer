FROM node:current-alpine
WORKDIR /app
COPY package.json pnpm-lock.yaml* ./
RUN npm install -g pnpm@9 && pnpm install --frozen-lockfile
COPY . .
RUN pnpm run build && pnpm prune --prod
ENV NODE_ENV=production PORT=3000 DATASETS_DIR=/data/datasets PUBLICATIONS_DIR=/data/publications PLOT_DIR=/data/plots
EXPOSE 3000
CMD ["node","build/index.js"]
