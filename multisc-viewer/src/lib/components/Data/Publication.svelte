<script lang="ts">
	import type { Publications } from '$lib/types/daemon';

	let { publication }: { publication: Publications[string] } = $props();

	let name = $derived(publication.name ?? publication._id);
</script>

<section class="mx-auto max-w-256">
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
			{publication.date}
			<!-- eslint-disable svelte/no-navigation-without-resolve -->
			<a
				class="anchor"
				href={publication.url}
				target="_blank"
				title={`${name} publication page`}
				rel="external noopener noreferrer"
			>
				{publication.identifier}
			</a>
			<!-- eslint-enable svelte/no-navigation-without-resolve -->
		</span>
	</section>
	<section>
		<h2 class="h6">Abstract</h2>
		<!-- eslint-disable-next-line svelte/no-at-html-tags -->
		<p>{@html publication.description}</p>
	</section>
</section>
