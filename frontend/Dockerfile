FROM node:20-slim AS base

# Install dependencies only when needed
FROM base AS deps
WORKDIR /application

# Copy package files
COPY package.json package-lock.json ./
RUN npm ci

# Rebuild the source code only when needed
FROM base AS builder
WORKDIR /application

COPY --from=deps /application/node_modules ./node_modules
COPY . .

# Generate API client from OpenAPI spec
RUN npm run generate:api

# Debug: Verify generated files exist
RUN ls -la api/ && echo "Generated API files confirmed"

# Set environment variable for build
ENV NEXT_PUBLIC_API_URL=http://backend:8000

# Build Next.js app (without turbopack for Docker compatibility)
RUN npm run build

# Production image
FROM base AS runner
WORKDIR /application

ENV NODE_ENV=production

# Create nextjs user
RUN addgroup --system --gid 1001 nodejs
RUN adduser --system --uid 1001 nextjs

# Copy built application
COPY --from=builder /application/public ./public
COPY --from=builder --chown=nextjs:nodejs /application/.next/standalone ./
COPY --from=builder --chown=nextjs:nodejs /application/.next/static ./.next/static

USER nextjs

EXPOSE 3000

ENV PORT=3000
ENV HOSTNAME="0.0.0.0"

CMD ["node", "server.js"]
