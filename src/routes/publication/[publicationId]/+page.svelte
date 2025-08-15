<script lang="ts">
	import { enhance } from '$app/forms';
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
</script>

<svelte:head>
	<title>{data.publication.title}</title>
</svelte:head>

{#await data.publication}
	<p>Loading publication...</p>
{:then publication}
	<section class="mx-auto max-w-256">
		<h2 class="text-center h2">{publication.title}</h2>
	</section>
	<section class="mx-auto max-w-256">
		<p>{publication.abstract}</p>
	</section>

	<hr class="hr" />

	<section>
		<h2 class="text-center h2">Datasets</h2>
		<form method="POST" action="?/plot" use:enhance>
			<div>
				<button type="submit" class="btn preset-filled">Plot</button>
			</div>
			<DatasetsTable datasets={publication.datasets} />
		</form>
	</section>
{/await}
