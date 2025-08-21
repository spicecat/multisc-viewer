<script lang="ts">
	import type { Publication } from '$lib/types/data';
	import { type PublicationData, publicationColumns } from '$lib/types/data-table';
	import DataTableSearch from './DataTableSearch.svelte';

	let { publications }: { publications: Publication[] } = $props();

	let data: PublicationData[] = $derived(
		publications.map((pub) => ({
			...pub,
			datasets: pub.datasets.map((ds) => ds.title),
			href: `/publication/${pub.id}`
		}))
	);
</script>

<DataTableSearch name="publications" {data} columns={publicationColumns} />
