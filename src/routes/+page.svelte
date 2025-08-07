<script lang="ts">
	import { enhance } from '$app/forms';
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import PublicationsTable from '$lib/components/DataTable/PublicationsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
</script>

<svelte:head>
	<title>Home</title>
</svelte:head>

<section>
	<h2 class="text-center h2">MultiSC-Viewer</h2>
	<p>
		MultiSC-Viewer is a web-based tool designed to compare and visualize gene expression in multiple
		single cell/nucleus dataset side by side for any gene across different brain regions, disease
		conditions, and species.
	</p>
</section>

<section>
	<h2 class="text-center h3">Features</h2>
	<ul class="list-inside list-disc">
		<li>
			<b>Multiple datasets visualization and comparison:</b> MultiSC-Viewer is the only tool current
			available that supports a side-by-side comparative view of more than two datasets simultaneously.
		</li>
		<li>
			<b>Interactive Data Visualization:</b> Utilize dynamic, interactive graphs for a real-time exploration
			of datasets, facilitating a deeper understanding of cellular diversity and expression similarities
			and differences across datasets.
		</li>
	</ul>
</section>

<hr class="hr" />

<section>
	<h2 class="text-center h2">Publications</h2>
	<PublicationsTable publications={data.publications} />
</section>

<hr class="hr" />

<section>
	<h2 class="text-center h2">Datasets</h2>
	<form method="POST" action="?/plot" use:enhance>
		<div>
			<button type="submit" class="btn preset-filled">Plot</button>
		</div>
		<DatasetsTable datasets={data.datasets} />
	</form>
</section>
