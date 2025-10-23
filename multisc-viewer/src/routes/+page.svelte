<script lang="ts">
	import PublicationsTable from '$lib/components/Data/PublicationsTable.svelte';
	import PlotsForm from '$lib/components/Plots/Form.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { datasets, publications } = $derived(data);
</script>

<svelte:head>
	<title>{data.meta.title}</title>
	<meta name="description" content={data.meta.description} />
</svelte:head>

<section class="mx-auto max-w-5xl">
	<section>
		<h2 class="text-center h2">MultiSC-Viewer</h2>
		<p>
			MultiSC-Viewer is a web application for visualizing and comparing gene expression in multiple
			single cell/nucleus datasets across different brain regions, disease conditions, and species.
		</p>
	</section>
	<section>
		<h2 class="text-center h3">Features</h2>
		<ul class="list-inside list-disc">
			<li>
				<b>Multiple datasets visualization and comparison</b>: currently the only available tool
				that supports a side-by-side comparative view of more than two datasets simultaneously
			</li>
			<li>
				<b>Interactive plots</b>: interactive plots controls for real-time exploration of datasets
			</li>
			<li>
				<b>Plot streaming</b>: renders plots as they are generated, allowing users to view results
				without waiting for all plots to complete
			</li>
			<li>
				<b>Load balancing and scalability</b>: manages datasets across multiple plotting daemons for
				parallelized plot generation.
			</li>
		</ul>
	</section>
</section>

<hr class="hr" />

<section class="mx-auto max-w-5xl">
	<h2 class="text-center h2">Publications</h2>
	<PublicationsTable {datasets} {publications} />
</section>

<hr class="hr" />

<section class="mx-auto">
	<h2 class="text-center h2">Datasets</h2>
	<PlotsForm {datasets} />
</section>
