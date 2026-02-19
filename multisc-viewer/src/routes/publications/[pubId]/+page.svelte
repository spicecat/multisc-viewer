<script lang="ts">
import Publication from "$lib/components/Data/Publication.svelte";
import PlotsForm from "$lib/components/Plots/Form.svelte";
import type { PageProps } from "./$types";

const { data }: PageProps = $props();
const { datasets, publication } = $derived(data);
</script>

<svelte:head>
	<title>{data.meta.title}</title>
	<meta name="description" content={data.meta.description} />
</svelte:head>

<Publication {publication} />

<hr class="hr" />

<section class="mx-auto">
	<h2 class="text-center h2">Datasets</h2>
	<PlotsForm {datasets} />
</section>
