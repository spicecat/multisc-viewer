<script lang="ts">
import type { Publication } from "$lib/types/daemon";

const { publication }: { publication: Publication } = $props();

const name = $derived(publication.name ?? publication.id);
</script>

<section class="mx-auto max-w-5xl">
	<section>
		<h1 class="h4">
			<!-- eslint-disable-next-line svelte/no-at-html-tags -->
			{@html name}
		</h1>
		<span class="block">
			{publication.author?.map((a) => a.name).join(', ')}
		</span>
		<span class="inline-block">
			{publication.journalName}
			{publication.publicationDate}
			<!-- eslint-disable svelte/no-navigation-without-resolve -->
			{#if publication.url}
				<a
					href={publication.url}
					class="anchor"
					target="_blank"
					title={`${name} publication source`}
					rel="external noopener noreferrer"
				>
					{publication.identifier ?? publication.url}
				</a>
			{/if}
			<!-- eslint-enable svelte/no-navigation-without-resolve -->
		</span>
	</section>
	<section>
		<h2 class="h6">Abstract</h2>
		<!-- eslint-disable-next-line svelte/no-at-html-tags -->
		<p>{@html publication.description}</p>
	</section>
</section>
