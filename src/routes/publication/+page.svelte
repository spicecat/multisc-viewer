<script lang="ts">
	import PublicationsTable from '$lib/components/DataTable/PublicationsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
</script>

<svelte:head>
	<title>Publications</title>
</svelte:head>

<section>
	<h1 class="text-center h1">Publications</h1>
	<PublicationsTable publications={data.publications} />
</section>
