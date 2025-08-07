<script lang="ts">
	import { enhance } from '$app/forms';
	import DatasetsTable from '$lib/components/DataTable/DatasetsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
</script>

<svelte:head>
	<title>{data.publication.title}</title>
</svelte:head>

<section>
	<h2 class="text-center h2">{data.publication.title}</h2>
</section>
<section>
	<p>{data.publication.abstract}</p>
</section>

<hr class="hr" />

<section>
	<h2 class="text-center h2">Datasets</h2>
	<form method="POST" action="?/plot" use:enhance>
		<div>
			<button type="submit" class="btn preset-filled">Plot</button>
		</div>
		<DatasetsTable datasets={data.publication.datasets} />
	</form>
</section>
