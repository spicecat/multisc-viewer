<script lang="ts">
	import PublicationsTable from '$lib/components/DataTable/PublicationsTable.svelte';
	import type { PageProps } from './$types';

	let { data }: PageProps = $props();
	let { publications } = data;
</script>

<svelte:head>
	<title>Publications</title>
</svelte:head>

<section class="mx-auto size-fit">
	<h2 class="text-center h2">Publications</h2>
	<PublicationsTable {publications} />
</section>
