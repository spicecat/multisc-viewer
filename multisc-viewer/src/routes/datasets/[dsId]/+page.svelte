<script lang="ts">
import OntologyTerm from "$lib/components/Data/OntologyTerm.svelte";
import type { PageProps } from "./$types";

const { data }: PageProps = $props();
const { dataset } = $derived(data);
</script>

<svelte:head>
	<title>{data.meta.title}</title>
	<meta name="description" content={data.meta.description} />
</svelte:head>

<section class="mx-auto max-w-5xl">
	<section>
		<h1 class="h4">
			<!-- eslint-disable-next-line svelte/no-at-html-tags -->
			{@html dataset.name}
		</h1>
		<span class="block">
			{dataset.author?.map((a) => a.name).join(', ')}
		</span>
		<span class="inline-block">
			{dataset.publicationDate}
			<!-- eslint-disable svelte/no-navigation-without-resolve -->
			{#if dataset.url}
				<a
					class="anchor"
					href={dataset.url}
					target="_blank"
					title={`${dataset.displayName ?? dataset.name ?? dataset.id} dataset source`}
					rel="external noopener noreferrer"
				>
					{dataset.identifier ?? dataset.url}
				</a>
			{/if}
			<!-- eslint-enable svelte/no-navigation-without-resolve -->
		</span>
	</section>
	<section>
		<ul class="list-inside space-y-1">
			<li>
				<OntologyTerm ontologyTerm={dataset.cellType}><b>Cell Type: </b></OntologyTerm>
			</li>
			<li>
				<OntologyTerm ontologyTerm={dataset.tissue}><b>Tissue: </b></OntologyTerm>
			</li>
			<li>
				<OntologyTerm ontologyTerm={dataset.healthCondition}><b>Health Condition: </b></OntologyTerm
				>
			</li>
			<li>
				<OntologyTerm ontologyTerm={dataset.species}><b>Species: </b></OntologyTerm>
			</li>
		</ul>
	</section>
	<section>
		<h2 class="h6">Abstract</h2>
		<!-- eslint-disable-next-line svelte/no-at-html-tags -->
		<p>{@html dataset.description}</p>
	</section>
	<section>
		{#if dataset.citation?.length}
			<h2 class="h6">Citations</h2>
			<ul class="list-inside list-disc">
				{#each dataset.citation as cit, i (`li-cit-${i}`)}
					{@const citName = cit.name ?? cit.identifier ?? cit.url}
					<li>
						{citName}
						<!-- eslint-disable svelte/no-navigation-without-resolve -->
						<a
							class="anchor"
							href={cit.url}
							target="_blank"
							title={`${citName} citation source`}
							rel="external noopener noreferrer"
						>
							{cit.identifier}
						</a>
						<!-- eslint-enable svelte/no-navigation-without-resolve -->
					</li>
				{/each}
			</ul>
		{/if}
	</section>
</section>
