<script lang="ts">
	/* eslint svelte/no-at-html-tags: "off" */
	import PlotForm from '$lib/components/PlotForm/PlotForm.svelte';
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
<section class="mx-auto size-fit">
	<h2 class="text-center h2">Datasets</h2>
	<PlotForm datasets={publication.datasets} />
</section>
