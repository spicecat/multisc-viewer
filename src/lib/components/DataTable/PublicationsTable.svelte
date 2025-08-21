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

<section class="mx-auto size-fit">
	<h2 class="text-center h2">Publications</h2>
	<DataTableSearch name="publications" {data} columns={publicationColumns} />
</section>
