<script lang="ts">
  import TopAppBar, { Row, Section, Title } from "@smui/top-app-bar";
  import Button from "@smui/button";
</script>

<TopAppBar
  variant="standard"
  dense
  style="top: 0; left: 0; position: relative; margin-bottom: 1rem;"
>
  <Row>
    <Section>
      <Button href="/">
        <Title>Plot Viewer</Title>
      </Button>
    </Section>
    <Section align="end">
      <Button href="/about">
        <Title>About</Title>
      </Button>
    </Section>
  </Row>
</TopAppBar>
