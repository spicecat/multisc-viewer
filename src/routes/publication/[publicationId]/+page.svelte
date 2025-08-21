<script lang="ts">
	/* eslint svelte/no-at-html-tags: "off" */
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { publication } = data;
</script>

<svelte:head>
	<title>{data.publication.title}</title>
</svelte:head>

<section class="mx-auto max-w-256">
	<h2 class="h2">{publication.title}</h2>
</section>
<section class="mx-auto max-w-256">
	<p>{publication.year}</p>
	<p class="font-bold">{publication.authors.join(', ')}</p>
	<p>PMID: {publication.PMID}</p>
</section>
<section class="mx-auto max-w-256">
	<h3 class="h3">Abstract</h3>
	<p>{@html publication.abstract}</p>
</section>

<hr class="hr" />
<DatasetsTable datasets={publication.datasets} publicationId={publication.id} />
