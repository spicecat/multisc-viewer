import type { Publication } from '../interfaces/types';

/**
 * Configuration for available publications
 */
export const publications: Publication[] = [
  {
    publicationId: 'astro-ad',
    name: 'Astrocyte Publications in Alzheimer\'s Disease',
    description: 'Collection of publications examining astrocyte populations in Alzheimer\'s Disease across different brain regions',
    datasets: []  // Populated dynamically from app.service.ts
  },
  {
    publicationId: 'astro-hd',
    name: 'Astrocyte Publications in Huntington\'s Disease',
    description: 'Collection of publications examining astrocyte populations in Huntington\'s Disease across different brain regions',
    datasets: []  // Populated dynamically from app.service.ts
  }
];